"""
Görsel kontrol: kayıtta kullanılan yörüngeyi, drone görünür hâlde ve pencereli olarak tekrar oynatır.

- Ana pencere: drone'u arkadan izleyen kamera (SpringArmChase).
- Alt pencereler (SubWindows): cam0 (sol), cam1 (sağ), cam0 derinliği.
- Ekranın sol üstünde canlı bilgi: zemine göre yükseklik, hız, kamera eğimi, taban çizgisi.
- Planlanan yol dünyaya kırmızı çizgi olarak çizilir.

Yörünge, kayit.py'nin `yorunge()` fonksiyonuyla ve kayıttaki kayit_meta.json argümanlarıyla
üretilir; yani gösterilen yol veri setini üreten yolun aynısıdır.

    python gosteri.py --settings_yaz            # Multirotor gösteri ayarlarını yazar
    python gosteri.py --meta <sekans>/kayit_meta.json --tur 2
"""
import argparse
import json
import math
import os
import sys
import time
from types import SimpleNamespace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from kayit import yorunge  # noqa: E402

AYAR = os.path.expanduser("~/Documents/AirSim/settings.json")


def ayarlar(meta):
    a = meta["arguman"]
    W, H, fov, B, pitch = a["genislik_px"], a["yukseklik_px"], a["fov"], a["taban"], a["pitch"]
    cap = [{"ImageType": 0, "Width": W, "Height": H, "FOV_Degrees": fov, "MotionBlurAmount": 0},
           {"ImageType": 3, "Width": 512, "Height": 384, "FOV_Degrees": fov}]
    # Gösteri: kameralar gövdenin altına/önüne alınır (X=+0.30, Z=+0.12 m, NED), yoksa görünür drone
    # gövdesi görüşü kapatır. Kayıtta (ComputerVision modu) gövde yoktur, kameralar X=Z=0'dadır.
    cam = lambda y: {"CaptureSettings": cap, "X": 0.30, "Y": y, "Z": 0.12,
                     "Pitch": pitch, "Roll": 0.0, "Yaw": 0.0}
    return {
        "SettingsVersion": 1.2, "SimMode": "Multirotor", "ViewMode": "SpringArmChase",
        "ClockSpeed": 1.0,
        "SubWindows": [
            {"WindowID": 0, "CameraName": "cam0", "ImageType": 0, "VehicleName": "Drone1", "Visible": True},
            {"WindowID": 1, "CameraName": "cam1", "ImageType": 0, "VehicleName": "Drone1", "Visible": True},
            {"WindowID": 2, "CameraName": "cam0", "ImageType": 3, "VehicleName": "Drone1", "Visible": True}],
        "Vehicles": {"Drone1": {"VehicleType": "SimpleFlight", "AutoCreate": True,
                                "Cameras": {"cam0": cam(-B / 2), "cam1": cam(B / 2)}}},
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta", default=os.path.expanduser("~/sim_envs/sdata_airsim/NH_egik_h5/kayit_meta.json"))
    ap.add_argument("--settings_yaz", action="store_true")
    ap.add_argument("--tur", type=int, default=1)
    ap.add_argument("--hiz_carpani", type=float, default=1.0, help="1 = gerçek zaman")
    ap.add_argument("--kip", default="ucus", choices=["ucus", "poz"],
                    help="ucus: otopilotla gerçek uçuş (akıcı); poz: kayıt pozlarını zorla oynatır (titrer)")
    a = ap.parse_args()
    meta = json.load(open(a.meta))
    if a.settings_yaz:
        json.dump(ayarlar(meta), open(AYAR, "w"), indent=2)
        print("yazıldı:", AYAR)
        return

    sys.path.insert(0, os.path.join(HERE, "lib_ms"))
    import airsim
    c = airsim.MultirotorClient()
    c.confirmConnection()
    ka = SimpleNamespace(**meta["arguman"])
    yol = yorunge(ka)
    dt = 1.0 / ka.hz
    B, pitch = ka.taban, ka.pitch

    c.simFlushPersistentMarkers()
    c.simPlotLineStrip([airsim.Vector3r(float(p[0]), float(p[1]), float(p[2])) for p, _ in yol[::10]],
                       color_rgba=[1, 0, 0, 1], thickness=8, is_persistent=True)
    if a.kip == "ucus":
        # GERÇEK UÇUŞ: simple_flight otopilotu, fizik motoru açık. Poz zorlanmaz -> titreme yok.
        c.enableApiControl(True, "Drone1"); c.armDisarm(True, "Drone1")
        c.takeoffAsync(vehicle_name="Drone1").join()
        p0 = yol[0][0]
        c.moveToPositionAsync(float(p0[0]), float(p0[1]), float(p0[2]), 2.0, vehicle_name="Drone1").join()
        wp = [airsim.Vector3r(float(p[0]), float(p[1]), float(p[2])) for p, _ in yol[::20]] + \
             [airsim.Vector3r(*map(float, yol[-1][0]))]
        for tur in range(a.tur):
            f = c.moveOnPathAsync(wp, ka.hiz, 3600, airsim.DrivetrainType.ForwardOnly,
                                  airsim.YawMode(False, 0), -1, 1, vehicle_name="Drone1")
            t0 = time.time()
            while time.time() - t0 < len(yol) * dt * 1.5 + 20:
                k = c.getMultirotorState("Drone1").kinematics_estimated
                v = math.hypot(k.linear_velocity.x_val, k.linear_velocity.y_val)
                col = c.simGetCollisionInfo("Drone1")
                c.simPrintLogMessage("Kip: ", f"GERCEK UCUS (otopilot)  tur {tur + 1}/{a.tur}", 0)
                c.simPrintLogMessage("Zemine gore yukseklik [m]: ", f"{-k.position.z_val:.2f}", 0)
                c.simPrintLogMessage("Yatay hiz [m/s]: ", f"{v:.2f}", 0)
                c.simPrintLogMessage("Kamera egimi (ufka gore): ", f"{pitch:.0f} derece asagi", 0)
                c.simPrintLogMessage("Stereo taban cizgisi [m]: ", f"{B:.2f}  (cam0 sol, cam1 sag)", 0)
                c.simPrintLogMessage("Carpisma: ", f"{col.has_collided} {col.object_name}", 0)
                if np.linalg.norm(np.array([k.position.x_val, k.position.y_val]) - yol[-1][0][:2]) < 1.0 \
                        and time.time() - t0 > 10:
                    break
                time.sleep(0.2)
            c.moveToPositionAsync(float(p0[0]), float(p0[1]), float(p0[2]), ka.hiz, vehicle_name="Drone1").join()
        c.hoverAsync("Drone1").join()
        print("bitti")
        return

    for tur in range(a.tur):
        t0 = time.time()
        for i, (p, q) in enumerate(yol):
            c.simSetVehiclePose(airsim.Pose(airsim.Vector3r(*map(float, p)),
                                            airsim.Quaternionr(float(q[1]), float(q[2]), float(q[3]), float(q[0]))),
                                True, "Drone1")
            if i % 10 == 0:
                if i + 10 < len(yol):
                    v = np.linalg.norm(yol[i + 10][0][:2] - p[:2]) / (10 * dt)
                c.simPrintLogMessage("Kare ", f"{i}/{len(yol)}  tur {tur + 1}/{a.tur}", 0)
                c.simPrintLogMessage("Zemine gore yukseklik [m]: ", f"{-p[2]:.2f}", 0)
                c.simPrintLogMessage("Yatay hiz [m/s]: ", f"{v:.2f}", 0)
                c.simPrintLogMessage("Kamera egimi (ufka gore): ", f"{pitch:.0f} derece asagi", 0)
                c.simPrintLogMessage("Stereo taban cizgisi [m]: ", f"{B:.2f}  (cam0 sol, cam1 sag)", 0)
            bekle = t0 + (i + 1) * dt / a.hiz_carpani - time.time()
            if bekle > 0:
                time.sleep(bekle)
    print("bitti")


if __name__ == "__main__":
    main()
