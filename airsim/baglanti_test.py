import sys; sys.path.insert(0,'lib_ms')
import airsim, time
c=airsim.VehicleClient(timeout_value=20); print('baglaniyor',flush=True); c.confirmConnection()
print('pose',c.simGetVehiclePose().position, flush=True)
t=time.time(); r=c.simGetImages([airsim.ImageRequest('cam0',airsim.ImageType.Scene,False,False)]); print('img',r[0].width,r[0].height,len(r[0].image_data_uint8), round(time.time()-t,2),flush=True)
