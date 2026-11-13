import sys, time; sys.path.insert(0,'lib_ms')
import airsim, numpy as np, cv2
c=airsim.VehicleClient(timeout_value=30); c.confirmConnection()
def grab():
    r=c.simGetImages([airsim.ImageRequest('cam0',airsim.ImageType.Scene,False,False)])[0]
    return cv2.cvtColor(np.frombuffer(r.image_data_uint8,np.uint8).reshape(r.height,r.width,-1)[:,:,:3],cv2.COLOR_BGR2GRAY).astype(float)
c.simPause(True)
c.simSetVehiclePose(airsim.Pose(airsim.Vector3r(8,-0.2,-5), airsim.to_quaternion(0,0,0)),True); time.sleep(0.5)
a=grab(); time.sleep(3); b=grab()
dd=np.abs(a-b); print('ort fark',dd.mean().round(2),' >20 piksel orani',(dd>20).mean().round(4),' >50',(dd>50).mean().round(4))
cv2.imwrite('<tmp>/fark.jpg', cv2.resize(np.hstack([a, np.clip(dd*5,0,255)]).astype(np.uint8),(1024,384)))
for cmd in ['r.PostProcessAAQuality 0','r.DefaultFeature.AntiAliasing 0','r.EyeAdaptationQuality 0']:
    c.simRunConsoleCommand(cmd)
time.sleep(1); a=grab(); time.sleep(3); b=grab(); dd=np.abs(a-b)
print('AA+goz uyumu kapali: ort fark',dd.mean().round(2),' >20',(dd>20).mean().round(4))
c.simPause(False)
