"""
Ortam keşfi (F rejimi): dar açılı, hafif aşağı bakan stereo ile kalkış, 8 yöne dönüp
görüntü ve derinlik alır. Uzun sıyıran görüş hattı olan yönü bulmak için.

    python sahne_kesif.py --settings_yaz --fov 15 --pitch -3
    python sahne_kesif.py --irtifa 5 --cikti <dizin>
"""
import argparse
import json
import math
import os
import sys
import time

import cv2
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
AYAR = os.path.expanduser("~/Documents/AirSim/settings.json")


def ayarlar(fov, pitch, B, W=1024, H=768):
    cap = [{"ImageType": 0, "Width": W, "Height": H, "FOV_Degrees": fov, "MotionBlurAmount": 0},
           {"ImageType": 2, "Width": 256, "Height": 192, "FOV_Degrees": fov},
           {"ImageType": 3, "Width": 512, "Height": 384, "FOV_Degrees": fov}]
    # Gösteri/keşif: kameralar gövdenin önüne/altına alınır (görünür drone gövdesi görüşü kapatmasın).
    cam = lambda y: {"CaptureSettings": cap, "X": 0.30, "Y": y, "Z": 0.12, "Pitch": pitch, "Roll": 0.0, "Yaw": 0.0}
    genis = {"CaptureSettings": [{"ImageType": 0, "Width": 640, "Height": 480, "FOV_Degrees": 90}],
             "X": 0.30, "Y": 0.0, "Z": 0.12, "Pitch": pitch, "Roll": 0.0, "Yaw": 0.0}
    return {"SettingsVersion": 1.2, "SimMode": "Multirotor", "ViewMode": "SpringArmChase", "ClockSpeed": 1.0,
            "SubWindows": [
                {"WindowID": 0, "CameraName": "cam0", "ImageType": 0, "VehicleName": "Drone1", "Visible": True},
                {"WindowID": 1, "CameraName": "cam1", "ImageType": 0, "VehicleName": "Drone1", "Visible": True},
                {"WindowID": 2, "CameraName": "genis", "ImageType": 0, "VehicleName": "Drone1", "Visible": True}],
            "Vehicles": {"Drone1": {"VehicleType": "SimpleFlight", "AutoCreate": True,
                                    "Cameras": {"cam0": cam(-B / 2), "cam1": cam(B / 2), "genis": genis}}}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--settings_yaz", action="store_true")
    ap.add_argument("--fov", type=float, default=15.0)
    ap.add_argument("--pitch", type=float, default=-3.0)
    ap.add_argument("--taban", type=float, default=0.40)
    ap.add_argument("--irtifa", type=float, default=5.0)
    ap.add_argument("--cikti", default=os.path.join(HERE, "kesif"))
    a = ap.parse_args()
    if a.settings_yaz:
        json.dump(ayarlar(a.fov, a.pitch, a.taban), open(AYAR, "w"), indent=2)
        print("yazıldı:", AYAR); return

    sys.path.insert(0, os.path.join(HERE, "lib_ms"))
    import airsim
    os.makedirs(a.cikti, exist_ok=True)
    c = airsim.MultirotorClient(); c.confirmConnection()
    c.enableApiControl(True, "Drone1"); c.armDisarm(True, "Drone1")
    c.takeoffAsync(vehicle_name="Drone1").join()
    p = c.getMultirotorState("Drone1").kinematics_estimated.position
    c.moveToZAsync(-a.irtifa, 1.5, vehicle_name="Drone1").join(); time.sleep(2)
    ozet = []
    for yaw in range(0, 360, 45):
        c.rotateToYawAsync(yaw, vehicle_name="Drone1").join(); time.sleep(1.5)
        r = c.simGetImages([airsim.ImageRequest("cam0", airsim.ImageType.Scene, False, False),
                            airsim.ImageRequest("cam0", airsim.ImageType.DepthPerspective, True, False),
                            airsim.ImageRequest("genis", airsim.ImageType.Scene, False, False)], "Drone1")
        img = np.frombuffer(r[0].image_data_uint8, np.uint8).reshape(r[0].height, r[0].width, 3)
        Z = np.array(r[1].image_data_float, np.float32).reshape(r[1].height, r[1].width)
        gen = np.frombuffer(r[2].image_data_uint8, np.uint8).reshape(r[2].height, r[2].width, 3)
        zemin = Z[Z < 5000]
        g = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
        n_kp = len(cv2.ORB_create(2000).detect(g, None))
        s = {"yaw": yaw, "derinlik_p10": float(np.percentile(zemin, 10)) if zemin.size else None,
             "derinlik_p50": float(np.median(zemin)) if zemin.size else None,
             "derinlik_p90": float(np.percentile(zemin, 90)) if zemin.size else None,
             "gok_orani": float(np.mean(Z >= 5000)), "orb_kp": n_kp,
             "kontrast": float(g.std())}
        ozet.append(s); print(json.dumps(s))
        cv2.imwrite(os.path.join(a.cikti, f"tele_yaw{yaw:03d}.jpg"), img)
        cv2.imwrite(os.path.join(a.cikti, f"genis_yaw{yaw:03d}.jpg"), gen)
    json.dump({"baslangic": [p.x_val, p.y_val, p.z_val], "irtifa": a.irtifa, "fov": a.fov,
               "pitch": a.pitch, "yonler": ozet}, open(os.path.join(a.cikti, "ozet.json"), "w"), indent=2)


if __name__ == "__main__":
    main()
