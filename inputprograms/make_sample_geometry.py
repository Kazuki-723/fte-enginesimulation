import numpy as np
import matplotlib.pyplot as plt

def make_circle(d, origin, N:int):  # d:穴の直径, origin:円の中心, N:点の数 Nを小さくすれば多角形も作れる．
    geometry = np.zeros((N,2))
    for i in range(N):
        geometry[i] = [d/2*np.cos(2*np.pi*i/N), d/2*np.sin(2*np.pi*i/N)]
    geometry += origin
    print(f"analytical lp = {np.pi*d}")
    print(f"analytical Ap = {np.pi*d**2/4}")
    return geometry

def make_gizagiza(d, D, n:int): # d:穴の直径, D:(d+ギザの高さ)/2, n:ギザの数
    geometry = np.zeros((2*n,2))
    for i in range(n):
        geometry[2*i] = [D/2*np.sin(2*np.pi*(2*i)/(2*n)), D/2*np.cos(2*np.pi*(2*i)/(2*n))]
        geometry[2*i+1] = [d/2*np.sin(2*np.pi*(2*i+1)/(2*n)), d/2*np.cos(2*np.pi*(2*i+1)/(2*n))]
    print(f"analytical lp = {np.sqrt((D/2)**2+(d/2)**2-2*D/2*d/2*np.cos(np.pi/n))*2*n}")
    print(f"analytical Ap = {1/2*D/2*d/2*np.sin(np.pi/n)*2*n}")
    return geometry

def make_gear(d, D, ratio, n:int):
    geometry = np.zeros((4*n,2))
    for i in range(0,4*n,4):
        center = 2*np.pi*(i+2)/(4*n)
        theta = [center-2*np.pi*2/(4*n)*ratio[1]/sum(ratio), center+2*np.pi*2/(4*n)*ratio[1]/sum(ratio)]
        geometry[i] = [D/2*np.sin(theta[0]), D/2*np.cos(theta[0])]
        geometry[i+1] = [d/2*np.sin(theta[0]), d/2*np.cos(theta[0])]
        geometry[i+2] = [d/2*np.sin(theta[1]), d/2*np.cos(theta[1])]
        geometry[i+3] = [D/2*np.sin(theta[1]), D/2*np.cos(theta[1])]
    print(f"analytical lp = {np.sum(np.linalg.norm(np.append(geometry[1:],geometry[0][np.newaxis,:],axis=0) - geometry,axis=1),axis=0)}")
    return geometry

def koch_snowflake(order, scale=10):    # matplotのチュートリアルから丸パクリ
    """
    Return two lists x, y of point coordinates of the Koch snowflake.

    Parameters
    ----------
    order : int
        The recursion depth.
    scale : float
        The extent of the snowflake (edge length of the base triangle).
    """
    def _koch_snowflake_complex(order):
        if order == 0:
            # initial triangle
            angles = np.array([0, 120, 240]) + 90
            return scale / np.sqrt(3) * np.exp(np.deg2rad(angles) * 1j)
        else:
            ZR = 0.5 - 0.5j * np.sqrt(3) / 3

            p1 = _koch_snowflake_complex(order - 1)  # start points
            p2 = np.roll(p1, shift=-1)  # end points
            dp = p2 - p1  # connection vectors

            new_points = np.empty(len(p1) * 4, dtype=np.complex128)
            new_points[::4] = p1
            new_points[1::4] = p1 + dp / 3
            new_points[2::4] = p1 + dp * ZR
            new_points[3::4] = p1 + dp / 3 * 2
            return new_points

    points = _koch_snowflake_complex(order)
    x, y = points.real, points.imag
    return x, y

def interpolate_geometry(geometry, N):  # geometryの点の間を直線で補完する．全体で点がN個になるように
    N -= len(geometry)
    M = N//len(geometry)
    if M < 1:
        print("N must be larger than length of geometry*2 ")
        exit(1)
    geometry_new = np.zeros(((M+1)*len(geometry)+1,2))
    geometry = np.append(geometry, geometry[0,np.newaxis], axis=0)
    for i in range(len(geometry)-1):
        geometry_new[(M+1)*i:(M+1)*(i+1)+1,0] = np.linspace(geometry[i,0],geometry[i+1,0],M+2)
        geometry_new[(M+1)*i:(M+1)*(i+1)+1,1] = np.linspace(geometry[i,1],geometry[i+1,1],M+2)
    geometry_new = np.delete(geometry_new,[-1],axis=0)
    print(f"len(geometry_new)={len(geometry_new)}")
    return geometry_new

def dh_maxmin(geometry):
    v_lines = np.append(geometry[1:],geometry[0][np.newaxis,:],axis=0) - geometry
    dh = np.linalg.norm(v_lines,axis=1)
    print(f"max : {np.max(dh)}")
    print(f"min : {np.min(dh)}")
    

if __name__=='__main__':
    #geometry = make_circle(0.034, [0,0], 360*1) # d[m], origin, 
    #geometry = make_gizagiza(0.02, 0.04, 4)
    geometry = make_gear(0.02, 0.04, (1,4), 8)
    geometry = interpolate_geometry(geometry, 360)
    # 点の間の距離の最大値最小値
    dh_maxmin(geometry)
    # 描画
    fig, ax = plt.subplots()
    ax.scatter(geometry[:,0], geometry[:,1], s=1, label="geometry")
    plt.show()
    # 保存
    np.savetxt("sample_geometry.csv", geometry, fmt='%.6f', delimiter=",")