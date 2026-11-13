import sys; sys.path.insert(0,'lib')
import airsim, numpy as np, cv2, time
c=airsim.MultirotorClient(); c.confirmConnection()
c.enableApiControl(True); c.armDisarm(True)
c.takeoffAsync().join(); c.moveToZAsync(-float(sys.argv[1]) if len(sys.argv)>1 else -15,3).join(); time.sleep(8)
r=c.simGetImages([airsim.ImageRequest('cam0',airsim.ImageType.Scene,False,False),
                  airsim.ImageRequest('cam1',airsim.ImageType.Scene,False,False),
                  airsim.ImageRequest('cam0',airsim.ImageType.DepthPlanar,True,False)])
for i,x in enumerate(r[:2]):
    img=np.frombuffer(x.image_data_uint8,np.uint8).reshape(x.height,x.width,-1)
    cv2.imwrite(f'duman/asagi_cam{i}.png',img); print('img',i,img.mean(),img.std())
d=airsim.list_to_2d_float_array(r[2].image_data_float,r[2].width,r[2].height)
print('derinlik p1/p50/p99',np.percentile(d,[1,50,99]))
