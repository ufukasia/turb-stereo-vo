import sys, time, math; sys.path.insert(0,'lib_ms')
import airsim, numpy as np
c=airsim.VehicleClient(timeout_value=30); c.confirmConnection()
def grab():
    r=c.simGetImages([airsim.ImageRequest('cam0',airsim.ImageType.Scene,False,False)])[0]
    return np.frombuffer(r.image_data_uint8,np.uint8).reshape(r.height,r.width,-1)[:,:,:3].astype(float)
P=lambda x: airsim.Pose(airsim.Vector3r(x,-0.2,-5), airsim.to_quaternion(0,0,0))
for mod in ['calisan','duraklatilmis']:
    c.simPause(mod=='duraklatilmis')
    c.simSetVehiclePose(P(8.0),True); time.sleep(0.5); a=grab(); time.sleep(2.0); b=grab()
    c.simSetVehiclePose(P(9.0),True); time.sleep(0.3); d=grab()
    print(mod,'ayni poz 2 s arayla fark',round(np.abs(a-b).mean(),3),' poz degisince fark',round(np.abs(a-d).mean(),2))
c.simPause(False)
