import sys, time; sys.path.insert(0,'lib_ms')
import airsim
c=airsim.VehicleClient(timeout_value=60); c.confirmConnection()
T=airsim.ImageType
for name,req in [('cam0 scene',[('cam0',T.Scene,False)]),('cam1 scene',[('cam1',T.Scene,False)]),('cam0+cam1 scene',[('cam0',T.Scene,False),('cam1',T.Scene,False)]),('cam0 depth',[('cam0',T.DepthPlanar,True)]),('cam0 scene x2',[('cam0',T.Scene,False)])]:
    t=time.time(); r=c.simGetImages([airsim.ImageRequest(a,b,f,False) for a,b,f in req]); print(name, round(time.time()-t,2), flush=True)
