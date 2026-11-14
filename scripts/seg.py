import cv2,numpy as np,math
def seg(path,S,lo=176,it=2):
    im=cv2.imread(path,0)
    c=im[311:1170,850:1714]
    up=cv2.resize(c,None,fx=3,fy=3,interpolation=cv2.INTER_CUBIC)
    fill=((up>=lo)&(up<=190)).astype(np.uint8)
    fill=cv2.erode(fill,np.ones((3,3),np.uint8),iterations=it)
    num,lab,st,cen=cv2.connectedComponentsWithStats(fill,connectivity=4)
    exp=(857*3)**2/S**2; sc=857/S; out=[]
    for i in range(1,num):
        if st[i,4]<0.5*exp: continue
        pts=cv2.findContours((lab==i).astype(np.uint8),cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_NONE)[0][0]
        (cx,cy),(w,h),ang=cv2.minAreaRect(pts)
        a=math.radians(ang%90); a=a-math.pi/2 if a>math.pi/4 else a
        out.append(((cx/3-3)/sc,(cy/3-3)/sc,a))
    return np.array(out)
