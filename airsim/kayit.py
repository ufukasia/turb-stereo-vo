"""
AirSim (ComputerVision modu) -> EuRoC ASL (mav0) stereo kayıt aracı.

NEDEN ComputerVision modu: fizik motoru yok; pozu biz veririz. Böylece
  - GT tanım gereği tam ve kamera çerçevesindedir,
  - cam0, cam1 ve derinlik AYNI render karesinden gelir (senkron),
  - kayıt deterministtir (aynı yörünge = aynı görüntüler).
TartanAir de bu yöntemle üretilmiştir.

KAMERA GEOMETRİSİ (TASARIM.md §4):
  aşağı bakan (nadir, Pitch=-90) stereo çift, taban çizgisi B gövde +y ekseninde.
  Nadir kamerada görüntü-x = gövde +y olduğundan cam1, cam0'a göre görüntü +x yönündedir
  (EuRoC/TartanAir ile aynı işaret: T_cn_cnm1[0,3] = -B).

ÇERÇEVELER:
  AirSim dünya = NED; AirSim kamera ekseni: x = optik eksen, y = sağ, z = aşağı.
  Görüntü çerçevesi (OpenCV): x sağ, y aşağı, z optik eksen.
      R_cam_img = [[0,0,1],[1,0,0],[0,1,0]]   (img z -> cam x, img x -> cam y, img y -> cam z)
  GT, cam0'ın GÖRÜNTÜ çerçevesinin dünya pozu olarak yazılır; sensor.yaml'da T_BS = I
  (gövde = rektifiye cam0) — turb/datasets/tartanair.py ile aynı sözleşme.

ZAMAN: 19 haneli ns damga (sözlüksel sıra = sayısal sıra; TARTANAIR_SIRALAMA_KUSURU dersi).

Kullanım:
  python kayit.py --out /media/.../sdata_airsim --seq NH_tarama_h15 \
      --senaryo tarama --irtifa 15 --hiz 3 --sure 90
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

import numpy as np
import cv2

HERE = os.path.dirname(os.path.abspath(__file__))

R_CAM_IMG = np.array([[0.0, 0.0, 1.0],
                      [1.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0]])
T0_NS = 1_700_000_000_000_000_000          # 19 hane


# ----------------------------------------------------------------------------
# Yörünge üretimi: drone benzeri, düzgün (C2) yol + ivmeden türetilen eğim
# ----------------------------------------------------------------------------
def _yol_noktalari(senaryo: str, uzunluk: float, genislik: float, seritler: int):
    """Yatay düzlemde (x kuzey, y doğu) yol noktaları; başlangıç (0,0)."""
    if senaryo == "tarama":                 # biçerdöver (lawnmower)
        pts = []
        for i in range(seritler):
            y = i * genislik / max(seritler - 1, 1)
            xs = (0.0, uzunluk) if i % 2 == 0 else (uzunluk, 0.0)
            pts += [(xs[0], y), (xs[1], y)]
        return np.array(pts)
    if senaryo == "orbit":
        r = genislik / 2
        a = np.linspace(0, 2 * np.pi, 73)
        return np.c_[r * np.cos(a) - r, r * np.sin(a)]
    if senaryo == "duz":
        return np.array([(0.0, 0.0), (uzunluk, 0.0)])
    raise ValueError(senaryo)


def _yeniden_ornekle(pts: np.ndarray, hiz: float, dt: float, kose_yaricapi: float):
    """Köşeleri yuvarlatılmış, sabit hızlı yol. Köşe yuvarlatma: Chaikin + yay uzunluğu."""
    p = pts.astype(float)
    for _ in range(4):                       # Chaikin köşe yumuşatma
        q = [p[0]]
        for a, b in zip(p[:-1], p[1:]):
            q += [0.75 * a + 0.25 * b, 0.25 * a + 0.75 * b]
        q.append(p[-1])
        p = np.array(q)
    seg = np.linalg.norm(np.diff(p, axis=0), axis=1)
    s = np.r_[0, np.cumsum(seg)]
    n = int(s[-1] / (hiz * dt))
    st = np.linspace(0, s[-1], n)
    x = np.interp(st, s, p[:, 0]); y = np.interp(st, s, p[:, 1])
    # hafif düzgünleştirme (ivme sürekliliği)
    k = max(3, int(kose_yaricapi / (hiz * dt)) | 1)
    w = np.hanning(k + 2)[1:-1]; w /= w.sum()
    xp = np.pad(x, k, mode="edge"); yp = np.pad(y, k, mode="edge")
    x = np.convolve(xp, w, mode="same")[k:-k]; y = np.convolve(yp, w, mode="same")[k:-k]
    return x, y


def _rota_oku(yol):
    """rota.json: {"yol": [[x,y],...], "zemin": {"x0","y0","res","z": [[...]]}}  (NED, m).
    zemin.z: çıplak arazi (DTM) NED z'si; satır = y ekseni değil, satır i -> x0 + i*res, sütun j -> y0 + j*res."""
    r = json.load(open(yol))
    return np.array(r["yol"], float), r.get("zemin")


def _zemin_z(zemin, x, y):
    """DTM'den çift doğrusal örnekleme."""
    Z = np.asarray(zemin["z"], float)
    fi = (np.asarray(x) - zemin["x0"]) / zemin["res"]
    fj = (np.asarray(y) - zemin["y0"]) / zemin["res"]
    i0 = np.clip(np.floor(fi).astype(int), 0, Z.shape[0] - 2)
    j0 = np.clip(np.floor(fj).astype(int), 0, Z.shape[1] - 2)
    di = np.clip(fi - i0, 0, 1); dj = np.clip(fj - j0, 0, 1)
    return (Z[i0, j0] * (1 - di) * (1 - dj) + Z[i0 + 1, j0] * di * (1 - dj)
            + Z[i0, j0 + 1] * (1 - di) * dj + Z[i0 + 1, j0 + 1] * di * dj)


def yorunge(a) -> list[tuple[np.ndarray, np.ndarray]]:
    """(konum_NED, kuaterniyon_wxyz) listesi — gövde pozu."""
    dt = 1.0 / a.hz
    zemin = None
    if a.senaryo == "dosya":
        pts, zemin = _rota_oku(a.rota)
    else:
        pts = _yol_noktalari(a.senaryo, a.uzunluk, a.genislik, a.seritler)
    x, y = _yeniden_ornekle(pts, a.hiz, dt, a.kose)
    n = min(len(x), int(a.sure * a.hz)) if a.sure > 0 else len(x)
    x, y = x[:n], y[:n]
    vx = np.gradient(x, dt); vy = np.gradient(y, dt)
    ax = np.gradient(vx, dt); ay = np.gradient(vy, dt)
    yaw = np.unwrap(np.arctan2(vy, vx)) if a.yaw_hareket else np.zeros(n)
    g = 9.81
    out = []
    rng = np.random.default_rng(a.tohum)
    # küçük irtifa salınımı (gerçek uçuşta irtifa sabit kalmaz)
    dz = a.irtifa_salinim * np.sin(2 * np.pi * np.arange(n) * dt / 11.0 + rng.uniform(0, 6.28))
    # Arazi takibi: zemine göre sabit yükseklik; zemin profili ~4 s'lik pencereyle yumuşatılır
    # (gerçek bir drone tümsekleri birebir izlemez).
    if zemin is not None:
        zg = _zemin_z(zemin, x, y)
        k = max(3, int(4.0 / dt) | 1)
        w = np.hanning(k + 2)[1:-1]; w /= w.sum()
        zg = np.convolve(np.pad(zg, k, mode="edge"), w, mode="same")[k:-k]
    else:
        zg = np.full(n, a.z_zemin)
    for i in range(n):
        # ivmeyi gövde çerçevesine çevir -> quadrotor eğimi (küçük açı): pitch ~ -a_ileri/g
        c, s_ = math.cos(yaw[i]), math.sin(yaw[i])
        a_f = c * ax[i] + s_ * ay[i]
        a_r = -s_ * ax[i] + c * ay[i]
        pitch = -math.atan2(a_f, g)
        roll = math.atan2(a_r, g)
        q = _euler_to_quat(roll, pitch, yaw[i])
        pos = np.array([x[i], y[i], zg[i] - (a.irtifa + dz[i])])
        out.append((pos, q))
    return out


def _euler_to_quat(roll, pitch, yaw):
    cr, sr = math.cos(roll / 2), math.sin(roll / 2)
    cp, sp = math.cos(pitch / 2), math.sin(pitch / 2)
    cy, sy = math.cos(yaw / 2), math.sin(yaw / 2)
    return np.array([cr * cp * cy + sr * sp * sy,
                     sr * cp * cy - cr * sp * sy,
                     cr * sp * cy + sr * cp * sy,
                     cr * cp * sy - sr * sp * cy])


def _quat_to_R(q):                          # wxyz
    w, x, y, z = q
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                     [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                     [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])


def _R_to_quat(R):
    """Dönme matrisi -> wxyz (Shepperd yöntemi, scipy gerektirmez)."""
    t = np.trace(R)
    if t > 0:
        s = math.sqrt(t + 1.0) * 2
        q = [0.25 * s, (R[2, 1] - R[1, 2]) / s, (R[0, 2] - R[2, 0]) / s, (R[1, 0] - R[0, 1]) / s]
    elif R[0, 0] > R[1, 1] and R[0, 0] > R[2, 2]:
        s = math.sqrt(1.0 + R[0, 0] - R[1, 1] - R[2, 2]) * 2
        q = [(R[2, 1] - R[1, 2]) / s, 0.25 * s, (R[0, 1] + R[1, 0]) / s, (R[0, 2] + R[2, 0]) / s]
    elif R[1, 1] > R[2, 2]:
        s = math.sqrt(1.0 + R[1, 1] - R[0, 0] - R[2, 2]) * 2
        q = [(R[0, 2] - R[2, 0]) / s, (R[0, 1] + R[1, 0]) / s, 0.25 * s, (R[1, 2] + R[2, 1]) / s]
    else:
        s = math.sqrt(1.0 + R[2, 2] - R[0, 0] - R[1, 1]) * 2
        q = [(R[1, 0] - R[0, 1]) / s, (R[0, 2] + R[2, 0]) / s, (R[1, 2] + R[2, 1]) / s, 0.25 * s]
    q = np.array(q)
    return q / np.linalg.norm(q) * (1 if q[0] >= 0 else -1)


# ----------------------------------------------------------------------------
def settings_json(a) -> dict:
    def cam(y):
        cs = [{"ImageType": 0, "Width": a.genislik_px, "Height": a.yukseklik_px,
               "FOV_Degrees": a.fov, "MotionBlurAmount": 0, "AutoExposureSpeed": 100},
              {"ImageType": 1, "Width": a.derinlik_px[0], "Height": a.derinlik_px[1],
               "FOV_Degrees": a.fov}]
        return {"CaptureSettings": cs, "X": 0.0, "Y": y, "Z": 0.0,
                "Pitch": a.pitch, "Roll": 0.0, "Yaw": 0.0}
    return {"SettingsVersion": 1.2, "SimMode": "ComputerVision", "ViewMode": "NoDisplay",
            "ClockSpeed": 1.0,
            "CameraDefaults": {"CaptureSettings": []},
            "Vehicles": {"CV": {"VehicleType": "ComputerVision", "AutoCreate": True,
                                "Cameras": {"cam0": cam(-a.taban / 2), "cam1": cam(a.taban / 2),
                                            # yalnız harita.py kullanır (dik aşağı, yükseklik haritası)
                                            "harita": {"CaptureSettings": [
                                                {"ImageType": 0, "Width": 512, "Height": 512, "FOV_Degrees": 90},
                                                {"ImageType": 1, "Width": 512, "Height": 512, "FOV_Degrees": 90}],
                                                "X": 0.0, "Y": 0.0, "Z": 0.0, "Pitch": -90.0, "Roll": 0.0, "Yaw": 0.0}}}}}


def _sensor_yaml(path, name, T_BS, fx, cx, cy, w, h, hz, T_cn_cnm1=None):
    lines = ["sensor_type: camera", f"comment: AirSim -> ASL ({name})", "T_BS:",
             "  cols: 4", "  rows: 4",
             "  data: [" + ", ".join(f"{v:.9f}" for v in np.asarray(T_BS).ravel()) + "]",
             f"rate_hz: {hz}", f"resolution: [{w}, {h}]", "camera_model: pinhole",
             f"intrinsics: [{fx:.6f}, {fx:.6f}, {cx:.6f}, {cy:.6f}]",
             "distortion_model: radial-tangential",
             "distortion_coefficients: [0.0, 0.0, 0.0, 0.0]"]
    if T_cn_cnm1 is not None:
        lines.append("T_cn_cnm1: [" + ", ".join(f"{v:.9f}" for v in np.asarray(T_cn_cnm1).ravel()) + "]")
    open(path, "w").write("\n".join(lines) + "\n")


def temizle(c, desen):
    """Hareketli/animasyonlu ya da bakışa bağlı (su) nesneleri sahneden siler; silinenleri döndürür.
    Ortam her açılışta nesneleri geri getirir, bu yüzden her kayıtta çağrılır."""
    import re
    silinen = []
    for _ in range(5):                       # silme bazen ilk denemede tutmuyor; tekrar dene
        hedef = [n for n in c.simListSceneObjects(".*") if re.search(desen, n)]
        if not hedef:
            break
        silinen += [n for n in hedef if c.simDestroyObject(n)]
        time.sleep(0.5)
    kalan = [n for n in c.simListSceneObjects(".*") if re.search(desen, n)]
    if kalan:
        raise RuntimeError(f"silinemeyen nesneler: {kalan}")
    print(f"temizlendi: {len(silinen)} nesne", flush=True)
    return silinen


def kaydet(a):
    sys.path.insert(0, os.path.join(HERE, a.istemci))
    import airsim
    c = airsim.VehicleClient()
    c.confirmConnection()
    if a.pozlama_bias is not None:
        try:
            c.simRunConsoleCommand(f"r.DefaultFeature.AutoExposure.Bias {a.pozlama_bias}")
        except Exception:
            pass

    silinen = temizle(c, a.temizle) if a.temizle else []
    # Render kaynaklı bulanıklıkları kapat: sinematik alan derinliği (dar FOV'da çok güçlü ve
    # DERİNLİĞE BAĞLI — türbülans bulanıklığıyla karışır) ve hareket bulanıklığı.
    render_komut = []
    if a.render_temiz:
        for k in ("r.DepthOfFieldQuality 0", "r.MotionBlurQuality 0", "r.DefaultFeature.MotionBlur 0"):
            render_komut.append([k, bool(c.simRunConsoleCommand(k))])
        print("render ayarı:", render_komut, flush=True)
    traj = yorunge(a)
    seq = os.path.join(a.out, a.seq)
    for sub in ("mav0/cam0/data", "mav0/cam1/data", "mav0/depth0/data",
                "mav0/state_groundtruth_estimate0"):
        os.makedirs(os.path.join(seq, sub), exist_ok=True)

    W, H = a.genislik_px, a.yukseklik_px
    fx = W / 2 / math.tan(math.radians(a.fov / 2))
    cx, cy = W / 2, H / 2
    dt_ns = int(round(1e9 / a.hz))
    req = [airsim.ImageRequest("cam0", airsim.ImageType.Scene, False, False),
           airsim.ImageRequest("cam1", airsim.ImageType.Scene, False, False),
           airsim.ImageRequest("cam0", airsim.ImageType.DepthPlanar, True, False)]

    rows0, rows1, gt, kalite = [], [], [], []
    # ısınma: pozu ver, birkaç kare bekle (doku akışı + otomatik pozlama otursun)
    p0, q0 = traj[0]
    c.simSetVehiclePose(airsim.Pose(airsim.Vector3r(*p0), airsim.Quaternionr(q0[1], q0[2], q0[3], q0[0])), True)
    time.sleep(a.isinma)
    # Sahneyi dondur: rüzgârla sallanan bitki örtüsü ve gölgeleri katı olmayan hareket üretir
    # (NH_egik_h5 ilk kaydında ~6 m'den sonra poz geri sıçramalarının nedeni). Duraklatılmış
    # dünyada poz değiştirip render almak, statik sahne varsayımını sağlar.
    if a.dondur:
        c.simPause(True)
        time.sleep(0.5)
    t_start = time.time()
    for i, (p, q) in enumerate(traj):
        c.simSetVehiclePose(airsim.Pose(airsim.Vector3r(*p),
                                        airsim.Quaternionr(q[1], q[2], q[3], q[0])), True)
        r = c.simGetImages(req)
        if len(r) != 3 or r[0].width == 0:
            raise RuntimeError(f"kare {i}: görüntü alınamadı")
        ts = T0_NS + i * dt_ns
        name = f"{ts}.png"
        imgs = []
        for k in (0, 1):
            im = np.frombuffer(r[k].image_data_uint8, np.uint8).reshape(r[k].height, r[k].width, -1)
            g = cv2.cvtColor(im[:, :, :3], cv2.COLOR_BGR2GRAY)
            cv2.imwrite(os.path.join(seq, f"mav0/cam{k}/data", name), g)
            imgs.append(g)
        d = np.array(r[2].image_data_float, dtype=np.float32).reshape(r[2].height, r[2].width)
        np.save(os.path.join(seq, "mav0/depth0/data", f"{ts}.npy"), d.astype(np.float16))
        rows0.append(f"{ts},{name}"); rows1.append(f"{ts},{name}")
        # GT: cam0 görüntü çerçevesinin dünya pozu (render edilen kameranın kendi pozu)
        cp, cq = r[0].camera_position, r[0].camera_orientation
        R_w_cam = _quat_to_R([cq.w_val, cq.x_val, cq.y_val, cq.z_val])
        qi = _R_to_quat(R_w_cam @ R_CAM_IMG)
        gt.append(f"{ts},{cp.x_val:.9f},{cp.y_val:.9f},{cp.z_val:.9f},"
                  f"{qi[0]:.9f},{qi[1]:.9f},{qi[2]:.9f},{qi[3]:.9f}")
        fin = d[np.isfinite(d) & (d < 1e4)]
        kalite.append([float(imgs[0].mean()), float(imgs[0].std()),
                       float(np.median(fin)) if fin.size else float("nan"),
                       float((d >= 1e4).mean())])
        if i % 100 == 0:
            el = time.time() - t_start
            print(f"  {i}/{len(traj)}  {el:.0f} s  parlaklik={kalite[-1][0]:.0f} "
                  f"derinlik_med={kalite[-1][2]:.1f} m  bos={kalite[-1][3]:.2f}", flush=True)

    if a.dondur:
        c.simPause(False)
    hdr = "#timestamp [ns],filename\n"
    for k, rows in ((0, rows0), (1, rows1)):
        open(os.path.join(seq, f"mav0/cam{k}/data.csv"), "w").write(hdr + "\n".join(rows) + "\n")
    open(os.path.join(seq, "mav0/state_groundtruth_estimate0/data.csv"), "w").write(
        "#timestamp,p_RS_R_x [m],p_RS_R_y [m],p_RS_R_z [m],q_RS_w [],q_RS_x [],q_RS_y [],q_RS_z []\n"
        + "\n".join(gt) + "\n")
    T_BS = np.eye(4)
    T_cn = np.eye(4); T_cn[0, 3] = -a.taban
    T_BS1 = np.eye(4); T_BS1[0, 3] = a.taban
    _sensor_yaml(os.path.join(seq, "mav0/cam0/sensor.yaml"), "cam0", T_BS, fx, cx, cy, W, H, a.hz)
    _sensor_yaml(os.path.join(seq, "mav0/cam1/sensor.yaml"), "cam1", T_BS1, fx, cx, cy, W, H, a.hz, T_cn)
    kal = np.array(kalite)
    meta = {"arguman": vars(a), "fx": fx, "n_kare": len(traj), "silinen_nesneler": silinen, "render_komut": render_komut,
            "parlaklik_ort": float(kal[:, 0].mean()), "kontrast_ort": float(kal[:, 1].mean()),
            "derinlik_medyan_m": float(np.nanmedian(kal[:, 2])),
            "bos_derinlik_orani": float(kal[:, 3].mean()),
            "sure_s": time.time() - t_start}
    json.dump(meta, open(os.path.join(seq, "kayit_meta.json"), "w"), indent=2, ensure_ascii=False)
    print(json.dumps({k: v for k, v in meta.items() if k != "arguman"}, indent=2))


def main():
    ap = argparse.ArgumentParser(description="AirSim -> EuRoC ASL stereo kayıt")
    ap.add_argument("--out", required=True)
    ap.add_argument("--seq", required=True)
    ap.add_argument("--senaryo", default="tarama", choices=["tarama", "orbit", "duz", "dosya"])
    ap.add_argument("--rota", default=None, help="senaryo=dosya: rota.json (yol + DTM)")
    ap.add_argument("--temizle", default=None,
                    help="silinecek nesne adları için regex (Africa: 'Water_Plane|Poacher|Crocodile|Hippopotamus|Rhinoceros')")
    ap.add_argument("--irtifa", type=float, default=15.0, help="zemin üstü [m]")
    ap.add_argument("--z_zemin", type=float, default=0.0, help="başlangıç noktasına göre zemin NED z [m]")
    ap.add_argument("--irtifa_salinim", type=float, default=0.5)
    ap.add_argument("--hiz", type=float, default=3.0)
    ap.add_argument("--hz", type=float, default=20.0)
    ap.add_argument("--sure", type=float, default=90.0, help="<=0: yolun tamamı")
    ap.add_argument("--uzunluk", type=float, default=80.0)
    ap.add_argument("--genislik", type=float, default=40.0)
    ap.add_argument("--seritler", type=int, default=3)
    ap.add_argument("--kose", type=float, default=4.0)
    ap.add_argument("--yaw_hareket", type=int, default=1)
    ap.add_argument("--tohum", type=int, default=0)
    ap.add_argument("--genislik_px", type=int, default=1024)
    ap.add_argument("--yukseklik_px", type=int, default=768)
    ap.add_argument("--fov", type=float, default=70.0)
    ap.add_argument("--taban", type=float, default=0.40)
    ap.add_argument("--pitch", type=float, default=-90.0)
    ap.add_argument("--derinlik_px", type=int, nargs=2, default=[256, 192],
                    help="derinlik çözünürlüğü; Python msgpack-rpc büyük float dizilerinde çok yavaş "
                         "(1024x768 ~39 s/kare), türbülans yol uzunluğu için 256x192 yeterli")
    ap.add_argument("--isinma", type=float, default=5.0)
    ap.add_argument("--pozlama_bias", type=float, default=None)
    ap.add_argument("--dondur", type=int, default=1, help="1: simPause ile sahneyi dondur (rüzgâr/gölge)")
    ap.add_argument("--render_temiz", type=int, default=1,
                    help="1: alan derinliği ve hareket bulanıklığını konsol komutuyla kapat")
    ap.add_argument("--istemci", default="lib_ms", help="lib_ms (Microsoft) ya da lib (Colosseum)")
    ap.add_argument("--settings_yaz", action="store_true",
                    help="yalnız ~/Documents/AirSim/settings.json yaz ve çık")
    a = ap.parse_args()
    if a.settings_yaz:
        p = os.path.expanduser("~/Documents/AirSim/settings.json")
        json.dump(settings_json(a), open(p, "w"), indent=2)
        print("yazıldı:", p)
        return
    kaydet(a)


if __name__ == "__main__":
    main()
