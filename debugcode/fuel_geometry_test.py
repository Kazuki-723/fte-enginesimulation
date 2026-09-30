import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
import os
import sys
import time

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from inputprograms.importjson import JsoncLoader

# ============================================================
# 1. 従来法（全数探索）による距離場の計算
# ============================================================

def compute_levelset_bruteforce(grid_points, lines, geometry, N_x, N_y, max_x, min_x, max_y, min_y):
    # 点距離
    v_AP = grid_points[:,None,:] - lines[:,0][None,:,:]
    d_points = np.min(np.linalg.norm(v_AP, axis=2), axis=1)
    mask_near_border = d_points < 0.001

    # 線分距離補正（従来コード）
    grid_points_nb = grid_points[mask_near_border]

    v_AP = grid_points_nb[:,None,:] - lines[:,0][None,:,:]
    v_BP = grid_points_nb[:,None,:] - lines[:,1][None,:,:]
    v_AB = lines[:,1][None,:,:] - lines[:,0][None,:,:]

    mask = (np.sum(v_AP*v_AB,axis=2)>0) & (np.sum(v_BP*v_AB,axis=2)<0)

    param_A = lines[:,0,1] - lines[:,1,1]
    param_B = lines[:,1,0] - lines[:,0,0]
    param_C = lines[:,0,0]*lines[:,1,1] - lines[:,1,0]*lines[:,0,1]

    d2 = (param_A[None,:]*grid_points_nb[:,0][:,None] +
          param_B[None,:]*grid_points_nb[:,1][:,None] +
          param_C[None,:])**2 / (param_A**2 + param_B**2)

    big = ((max_x-min_x)*2 + (max_y-min_y)*2)
    d2 = d2*mask + big*(~mask)

    d_lines = np.sqrt(np.min(d2, axis=1))

    levelset = d_points.copy()
    levelset[mask_near_border] = np.minimum(d_points[mask_near_border], d_lines)
    levelset = levelset.reshape((N_x,N_y))
    
    polygon = Path(geometry)
    is_inside = polygon.contains_points(grid_points).reshape((N_x,N_y))
    levelset[is_inside]*=-1
    update = levelset < d_points.reshape((N_x,N_y))
    levelset = d_points.reshape((N_x,N_y))*~update + levelset*update

    return levelset


# ============================================================
# 2. KD-tree 最適化版による距離場の計算
# ============================================================

from scipy.spatial import cKDTree

def compute_levelset_kdtree(grid_points, lines, geometry, N_x, N_y, max_x, min_x, max_y, min_y):
    tree = cKDTree(geometry)

    # 点距離
    d_points, idx = tree.query(grid_points)
    mask_near_border = d_points < 0.001

    # 線分距離補正
    grid_points_nb = grid_points[mask_near_border]

    seg_idx1 = idx[mask_near_border]
    seg_idx2 = (seg_idx1 - 1) % len(lines)

    A = np.stack((lines[seg_idx1,0], lines[seg_idx2,0]), axis=1)
    B = np.stack((lines[seg_idx1,1], lines[seg_idx2,1]), axis=1)

    v_AP = grid_points_nb[:,None,:] - A
    v_BP = grid_points_nb[:,None,:] - B
    v_AB = B - A

    mask = (np.sum(v_AP*v_AB,axis=2)>0) & (np.sum(v_BP*v_AB,axis=2)<0)

    param_A = A[:,:,1] - B[:,:,1]
    param_B = B[:,:,0] - A[:,:,0]
    param_C = A[:,:,0]*B[:,:,1] - B[:,:,0]*A[:,:,1]

    d2 = (param_A*grid_points_nb[:,0][:,None] +
          param_B*grid_points_nb[:,1][:,None] +
          param_C)**2 / (param_A**2 + param_B**2)

    big = ((max_x-min_x)*2 + (max_y-min_y)*2)
    d2 = d2*mask + big*(~mask)

    d_lines = np.sqrt(np.min(d2, axis=1))

    levelset = d_points.copy()
    levelset[mask_near_border] = np.minimum(d_points[mask_near_border], d_lines)
    levelset = levelset.reshape((N_x,N_y))

    polygon = Path(geometry)
    is_inside = polygon.contains_points(grid_points).reshape((N_x,N_y))
    levelset[is_inside]*=-1
    update = levelset < d_points.reshape((N_x,N_y))
    levelset = d_points.reshape((N_x,N_y))*~update + levelset*update
   
    return levelset


# ============================================================
# 3. テスト実行
# ============================================================

# データ読み込み
# テスト対象のfilenameを入れる
setting_filename = "geometry_settings.jsonc"
try:
    loader = JsoncLoader(setting_filename)
    settings = loader.load()
except Exception as e:
    print(f"loading error: {e}")
    exit(1)

# 境界の読み込み
if settings["mode"]=="levelset":
    levelset = np.loadtxt(settings["levelset"]["filename"], delimiter=",", skiprows=1, dtype=float)
    N_x = len(levelset)
    N_y = len(levelset[0])
    min_x = settings["levelset"]["min_x"]
    max_x = settings["levelset"]["max_x"]
    min_y = settings["levelset"]["min_y"]
    max_y = settings["levelset"]["max_y"]
    delta_x = (max_x - min_x)/N_x
    delta_y = (max_y - min_y)/N_y
    symmetry = settings["levelset"]["symmetry"]

if settings["mode"]=="geometry":
    geometry_files = settings["geometry"]["geometry"]

    # 計算領域の作成
    min_x = settings["geometry"]["culc_area"]["min_x"]
    max_x = settings["geometry"]["culc_area"]["max_x"]
    N_x = settings["geometry"]["culc_area"]["N_x"]
    min_y = settings["geometry"]["culc_area"]["min_y"]
    max_y = settings["geometry"]["culc_area"]["max_y"]
    N_y = settings["geometry"]["culc_area"]["N_y"]
    delta_x = (max_x - min_x)/N_x
    delta_y = (max_y - min_y)/N_y
    levelset = np.zeros((N_x,N_y)) + (max_x-min_x)*2+(max_y-min_y)*2
    symmetry = 100

    for geometry_file in geometry_files:    # 穴が別ならファイルを分けているという仮定で．不便な気もする．一つのファイルで識別番号振らせるのとどっちがいいか．どっちもするべきか．
        # print(f"culclate {geometry_file["filename"]}")
        geometry = np.loadtxt(geometry_file["filename"], delimiter=',', dtype = float, encoding='utf-8')

        # この境界での計算領域のパラメータを作る．対称性を利用する計算のため．
        min_x = min_x
        max_x = max_x
        N_x = N_x
        min_y = min_y
        max_y = max_y
        N_y = N_y
        # 対称性から計算領域を狭める
        if geometry_file["symmetry"]==4:
            max_x = (max_x + min_x)/2
            N_x = int(N_x/2)
            max_y = (max_y + min_y)/2
            N_y = int(N_y/2)
        # 全体の対称性
        symmetry = min(symmetry, geometry_file["symmetry"])

        #phi = 0の生成
        lines = np.array([geometry, np.append(geometry[1:],geometry[0]).reshape(len(geometry),2)]).transpose(1,0,2)

        # Mesh生成
        x = np.linspace(min_x, max_x, N_x)
        y = np.linspace(min_y, max_y, N_y)
        X, Y = np.meshgrid(x, y)
        grid_points = np.column_stack((X.ravel(), Y.ravel()))

# テスト本体

# 従来型全探索
print("=== Running brute-force method ===")
t0 = time.perf_counter()
levelset_brute = compute_levelset_bruteforce(
    grid_points, lines, geometry, N_x, N_y, max_x, min_x, max_y, min_y)
t1 = time.perf_counter()
print(f"Brute-force time: {t1 - t0:.6f} sec")

# KD-treeによる探索
print("\n=== Running KD-tree method ===")
t2 = time.perf_counter()
levelset_kdtree = compute_levelset_kdtree(
    grid_points, lines, geometry, N_x, N_y, max_x, min_x, max_y, min_y)
t3 = time.perf_counter()
print(f"KD-tree time: {t3 - t2:.6f} sec")


# ============================================================
# 4. 差分の統計
# ============================================================

diff = levelset_brute - levelset_kdtree

print("=== KD-tree vs brute-force 差分統計 ===")
print("max abs diff:", np.max(np.abs(diff)))
print("mean abs diff:", np.mean(np.abs(diff)))
print("RMS diff:", np.sqrt(np.mean(diff**2)))

# ============================================================
# 5. 差分のヒストグラム
# ============================================================

plt.figure(figsize=(6,4))
plt.hist(diff.flatten(), bins=50)
plt.title("Difference Histogram (brute - KD-tree)")
plt.xlabel("difference")
plt.ylabel("count")
plt.tight_layout()
plt.show()

# ============================================================
# 6. 差分の空間分布
# ============================================================

plt.figure(figsize=(6,5))
plt.imshow(diff, origin='lower')
plt.colorbar()
plt.title("Spatial Distribution of Difference")
plt.tight_layout()
plt.show()

# ============================================================
# 7. それぞれの実行結果表示
# ============================================================
fig, ax = plt.subplots()
im  = ax.imshow(levelset_brute, vmin=np.min(levelset_brute), vmax=np.max(levelset_brute), cmap = "coolwarm")
cbar = fig.colorbar(im)
cbar.set_label("Distance From phi = 0", fontsize=10)
levels = np.arange(0,0.01,2e-3)
ctr = ax.contour(levelset_brute, levels, colors="black")
ax.clabel(ctr, levels, inline=1)
plt.title("Brute-force Levelset")
plt.tight_layout()

fig, ax = plt.subplots()
im  = ax.imshow(levelset_kdtree, vmin=np.min(levelset_kdtree), vmax=np.max(levelset_kdtree), cmap = "coolwarm")
cbar = fig.colorbar(im)
cbar.set_label("Distance From phi = 0", fontsize=10)
levels = np.arange(0,0.01,2e-3)
ctr = ax.contour(levelset_kdtree, levels, colors="black")
ax.clabel(ctr, levels, inline=1)
plt.title("KD-tree Levelset")
plt.tight_layout()
plt.show()