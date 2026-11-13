"""rota_sec.py adaylarını simülatörde değerlendirir: her şeridin başında, uçuş yönüne bakan
cam0'ın derinlik dağılımı (uzun sıyıran görüş hattı ölçütü). ComputerVision modu gerekir."""
import json, math, os, sys
from types import SimpleNamespace
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import kayit
from kayit import yorunge, _yol_noktalari, temizle
sys.path.insert(0, os.path.join(HERE, "lib_ms")); import airsim

def main(harita, rota_json, irtifa, cikti):
    d = np.load(harita); r = json.load(open(rota_json)); zemin = r["zemin"]; P = r["parametre"]
    c = airsim.VehicleClient(); c.confirmConnection()
    temizle(c, "Water_Plane|Poacher|Crocodile|Hippopotamus|Rhinoceros")
    gor, secilen = [], set()
    for a_ in r["adaylar"]:
        key = (round(a_["x"] / 40), round(a_["y"] / 40), a_["th"] // 20)
        if key in secilen: continue
        secilen.add(key)
        th = math.radians(a_["th"]); base = _yol_noktalari("tarama", P["L"], P["G"], 2)
        cs, sn = math.cos(th), math.sin(th)
        pts = np.c_[a_["x"] + cs * base[:, 0] - sn * base[:, 1], a_["y"] + sn * base[:, 0] + cs * base[:, 1]]
        kayit._rota_oku = lambda _: (pts, zemin)
        ka = SimpleNamespace(senaryo="dosya", rota=None, hz=20.0, hiz=P["hiz"], kose=P["kose"], sure=0.0, yaw_hareket=1,
                             tohum=0, irtifa=irtifa, irtifa_salinim=0.3, z_zemin=0.0, uzunluk=P["L"], genislik=P["G"], seritler=2)
        traj = yorunge(ka); n = len(traj)
        med = []
        for i in [int(0.1 * n), int(0.3 * n), int(0.7 * n), int(0.9 * n)]:
            p, q = traj[i]
            c.simSetVehiclePose(airsim.Pose(airsim.Vector3r(*map(float, p)), airsim.Quaternionr(float(q[1]), float(q[2]), float(q[3]), float(q[0]))), True)
            import time; time.sleep(0.6)
            rr = c.simGetImages([airsim.ImageRequest("cam0", airsim.ImageType.DepthPlanar, True, False)])[0]
            Z = np.array(rr.image_data_float, np.float32); Z = Z[Z < 5000]
            med.append(float(np.median(Z)) if Z.size else 0.0)
        s = dict(a_, derinlik_medyan=med, skor=min(med))
        gor.append(s); print(json.dumps(s), flush=True)
        if len(gor) >= 15: break
    gor.sort(key=lambda s: -s["skor"])
    json.dump(gor, open(cikti, "w"), indent=1)

if __name__ == "__main__":
    main(*sys.argv[1:3], float(sys.argv[3]), sys.argv[4])
