import sys, time; sys.path.insert(0,'lib_ms')
import airsim
c=airsim.VehicleClient(timeout_value=60); c.confirmConnection()
T=airsim.ImageType
for name,(cam,it,fl) in [('depth2',('cam0',T.DepthPlanar,True)),('depth3',('cam0',T.DepthPlanar,True)),('persp float',('cam0',T.DepthPerspective,True)),('depthvis png',('cam0',T.DepthVis,False)),('planar uint8',('cam0',T.DepthPlanar,False))]:
    t=time.time(); r=c.simGetImages([airsim.ImageRequest(cam,it,fl,False)]); print(name, round(time.time()-t,2), r[0].width, len(r[0].image_data_float), len(r[0].image_data_uint8), flush=True)
