import sys, math, time; sys.path.insert(0,'lib_ms')
import airsim, numpy as np, cv2
c=airsim.VehicleClient(timeout_value=30); c.confirmConnection()
tag0=sys.argv[1]; pts=[tuple(map(float,p.split(','))) for p in sys.argv[2:]]
for (x,y,z,yaw) in pts:
    c.simSetVehiclePose(airsim.Pose(airsim.Vector3r(x,y,z), airsim.to_quaternion(0,0,math.radians(yaw))),True); time.sleep(1.0)
    r=c.simGetImages([airsim.ImageRequest('cam0',airsim.ImageType.Scene,False,False),airsim.ImageRequest('cam0',airsim.ImageType.DepthPlanar,True,False)])
    im=np.frombuffer(r[0].image_data_uint8,np.uint8).reshape(r[0].height,r[0].width,-1)[:,:,:3]
    d=np.array(r[1].image_data_float,np.float32).reshape(r[1].height,r[1].width)
    tag=f'{tag0}_x{x:g}_y{y:g}_z{z:g}_yaw{yaw:g}'
    cv2.imwrite(f'duman/{tag}.jpg',cv2.resize(im,(512,384)))
    print(tag,'parl',round(im.mean(),1),'kontr',round(im.std(),1),'derinlik(alt/orta/ust)',[round(float(np.median(d[rr])),1) for rr in (int(d.shape[0]*0.96),d.shape[0]//2,int(d.shape[0]*0.08))],'bos',round(float((d>1e4).mean()),2),flush=True)
