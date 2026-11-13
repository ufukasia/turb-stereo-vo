import sys, math, time; sys.path.insert(0,'lib_ms')
import airsim, numpy as np
c=airsim.VehicleClient(timeout_value=15); c.confirmConnection()
for i,yaw in enumerate([180,270,0,45,135]):
    t=time.time(); c.simSetVehiclePose(airsim.Pose(airsim.Vector3r(0,0,-5), airsim.to_quaternion(0,0,math.radians(yaw))),True); t1=time.time()
    r=c.simGetImages([airsim.ImageRequest('cam0',airsim.ImageType.Scene,False,False),airsim.ImageRequest('cam1',airsim.ImageType.Scene,False,False),airsim.ImageRequest('cam0',airsim.ImageType.DepthPlanar,True,False)]); t2=time.time()
    print(i,yaw,'setpose',round(t1-t,2),'imgs',round(t2-t1,2),[x.width for x in r],flush=True)
