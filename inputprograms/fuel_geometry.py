import numpy as np
from scipy.ndimage import distance_transform_edt
from skimage.measure import find_contours
import matplotlib.pyplot as plt
from matplotlib.path import Path
import time
#from inputprograms.importjson import JsoncLoader
from importjson import JsoncLoader

# delta_x = delta_y を前提にしている．

class FuelGeometry:
    def __init__(self):
        self.r_arr=np.zeros(1)

    # 燃料端面phi = 0の周回長さ計算
    def culc_lp(self, levelset, symmetry=1):
        #if symmetry==4:
        #    levelset = levelset[int(len(levelset)/2):,int(len(levelset)/2):]
        ctr = find_contours(levelset, 0.0)[0] # phi=0の等高線
        lp = np.sum(np.linalg.norm(np.append(ctr[1:],ctr[0][np.newaxis,:],axis=0) - ctr,axis=1),axis=0)*self.delta_x
        #if symmetry==4:
        #    lp*=4
        return lp

    # 燃料端面phi = 0の内部の面積計算
    def culc_Ap(self, levelset, symmetry=1):
        # eps = self.delta_x*c
        #if self.symmetry==4:
        #    levelset = levelset[int(len(levelset)/2):,int(len(levelset)/2):]
        mask = levelset < 0
        A_p = np.sum(mask) * (self.delta_x*self.delta_y) #+ eps/2 + np.sum((abs(levelset) < eps)*(np.sin(levelset*np.pi/eps) + levelset))/2) * (self.delta_x*self.delta_y)
        #if self.symmetry==4:
        #    A_p*=4
        return A_p
    
    def culc_initial_levelset(self, setting_filename):
        # settingsの読み込み
        # print("input geometry settings file name:")
        self.setting_filename = setting_filename    # input("> ").strip()
        try:
            loader = JsoncLoader(self.setting_filename)
            settings = loader.load()
        except Exception as e:
            print(f"loading error: {e}")
            exit(1)

        # 実行時間計測
        start = time.perf_counter()
        # 境界の読み込み
        if settings["mode"]=="levelset":
            levelset = np.loadtxt(settings["levelset"]["filename"], delimiter=",", skiprows=1, dtype=float)
            self.N_x = len(levelset)
            self.N_y = len(levelset[0])
            self.min_x = settings["levelset"]["min_x"]
            self.max_x = settings["levelset"]["max_x"]
            self.min_y = settings["levelset"]["min_y"]
            self.max_y = settings["levelset"]["max_y"]
            self.delta_x = (self.max_x - self.min_x)/self.N_x
            self.delta_y = (self.max_y - self.min_y)/self.N_y
            self.symmetry = settings["levelset"]["symmetry"]

        if settings["mode"]=="geometry":
            geometry_files = settings["geometry"]["geometry"]
        
            # 計算領域の作成
            self.min_x = settings["geometry"]["culc_area"]["min_x"]
            self.max_x = settings["geometry"]["culc_area"]["max_x"]
            self.N_x = settings["geometry"]["culc_area"]["N_x"]
            self.min_y = settings["geometry"]["culc_area"]["min_y"]
            self.max_y = settings["geometry"]["culc_area"]["max_y"]
            self.N_y = settings["geometry"]["culc_area"]["N_y"]
            self.delta_x = (self.max_x - self.min_x)/self.N_x
            self.delta_y = (self.max_y - self.min_y)/self.N_y
            levelset = np.zeros((self.N_x,self.N_y)) + (self.max_x-self.min_x)*2+(self.max_y-self.min_y)*2
            self.symmetry = 100

            for geometry_file in geometry_files:    # 穴が別ならファイルを分けているという仮定で．不便な気もする．一つのファイルで識別番号振らせるのとどっちがいいか．どっちもするべきか．
                # print(f"culclate {geometry_file["filename"]}")
                geometry = np.loadtxt(geometry_file["filename"], delimiter=',', dtype = float, encoding='utf-8')

                # この境界での計算領域のパラメータを作る．対称性を利用する計算のため．
                min_x = self.min_x
                max_x = self.max_x
                N_x = self.N_x
                min_y = self.min_y
                max_y = self.max_y
                N_y = self.N_y
                # 対称性から計算領域を狭める
                if geometry_file["symmetry"]==4:
                    max_x = (self.max_x + self.min_x)/2
                    N_x = int(N_x/2)
                    max_y = (self.max_y + self.min_y)/2
                    N_y = int(N_y/2)
                # 全体の対称性
                self.symmetry = min(self.symmetry, geometry_file["symmetry"])

                #---------------------
                # levelset関数の計算
                #---------------------
                # Mesh生成
                x = np.linspace(min_x, max_x, N_x)
                y = np.linspace(min_y, max_y, N_y)
                X,Y = np.meshgrid(x,y)
                grid_points = np.column_stack((X.ravel(), Y.ravel()))
                # 境界線の読み込み，境界線上の点座標を保有
                lines = np.array([geometry, np.append(geometry[1:],geometry[0]).reshape(len(geometry),2)]).transpose(1,0,2)  # M行2列で各要素は1行2列(M,2,2)
                # linesの各点とgrid_pointsの各座標の差分(x,y)を計算する
                v_AP = grid_points[:,np.newaxis,:] - lines[:,0][np.newaxis,:,:]   # (N^2,1,2)+(1,M,2)->(N^2,M,2)
                # 格子点と点の距離
                # 全点計算
                d_points = np.linalg.norm(v_AP, axis=2)
                # 各gridに対する最小値の計算
                d_points = np.min(d_points, axis=1) #(N^2,1)

                # 距離関数の値をプロット
                fig, ax = plt.subplots()
                im  = ax.imshow(d_points.reshape((N_x,N_y)), vmin=min(d_points), vmax=max(d_points))
                cbar = fig.colorbar(im)
                cbar.set_label("Distance From phi = 0", fontsize=10)
                plt.show()

                # 境界の近くのみ線分との距離も計算
                mask_near_border = d_points < 0.001    # (N',1)

                # 距離関数が一定以下(上のmask)のみハイライトプロット
                plt.figure("mask_near_border")
                im  = plt.imshow(mask_near_border.reshape((N_x,N_y)))
                plt.show()

                # Meshの近傍点より近い点を探査する
                # 読み切れてないのでまた後で
                grid_points_nb = grid_points[mask_near_border]  #(N',1)
                v_AP = grid_points_nb[:,np.newaxis,:] - lines[:,0][np.newaxis,:,:]   # (N',1,2)+(1,M,2)->(N',M,2)
                v_BP = grid_points_nb[:,np.newaxis,:] - lines[:,1][np.newaxis,:,:]   # (N',M,2)
                v_AB = lines[:,1][np.newaxis,:,:] - lines[:,0][np.newaxis,:,:]  # (1,M,2)            
                # 線分の両端から格子点への角度がどちらも90°以下，つまり線分の両端より線分の方が近い場合を抽出
                mask = (np.sum(v_AP*v_AB,axis=2)>0)&(np.sum(v_BP*v_AB,axis=2)<0).astype(bool)  # (N',M)
                # 格子点と直線の距離 d^2 = (A*x+B*y+C)**2/(A**2+B**2)   A,B,C:linesから　x,y:pointsから
                param_lines = np.zeros((len(lines),3))  # (M,3)
                param_lines[:,0] = lines[:,0,1] - lines[:,1,1] # A=y1-y2
                param_lines[:,1] = lines[:,1,0] - lines[:,0,0] # B=x2-x1
                param_lines[:,2] = lines[:,0,0]*lines[:,1,1] - lines[:,1,0]*lines[:,0,1] # C=x1*y2-x2*y1
                d_lines = (param_lines[:,0]*grid_points_nb[:,0][:,np.newaxis]+param_lines[:,1]*grid_points_nb[:,1][:,np.newaxis]+param_lines[:,2])**2/ \
                (param_lines[:,0]**2+param_lines[:,1]**2)*mask + ((max_x-min_x)*2+(max_y-min_y)*2)*~mask # (N', M)
                d_lines = np.sqrt(np.min(d_lines, axis=1))   #(N',1)
                levelset_new = d_points.copy()
                levelset_new[mask_near_border] = np.min(np.stack((d_points[mask_near_border], d_lines), axis=1), axis=1)
                levelset_new = levelset_new.reshape((N_x,N_y))
                plt.figure("updated")
                im  = plt.imshow(levelset_new-d_points.reshape((N_x,N_y)))
                plt.show()
                # 符号付の値に変換，内部が負
                polygon = Path(geometry)
                is_inside = polygon.contains_points(grid_points).reshape((N_x,N_y))
                levelset_new[is_inside]*=-1
                # 対称性
                if geometry_file["symmetry"]==4:
                    levelset_new = np.concatenate((levelset_new, np.fliplr(levelset_new)), 1)
                    levelset_new = np.concatenate((levelset_new, np.flipud(levelset_new)), 0)
                # 更新
                update = levelset_new < levelset
                levelset = levelset*~update + levelset_new*update
        # 初期ポート断面積，周長の計算
        A_p_init = self.culc_Ap(levelset) # ポート断面積計算
        l_p_init = self.culc_lp(levelset) # 周回長さ計算

        # 実行時間 
        end = time.perf_counter()
        print(f"Elapsed time: {end - start:.6f} seconds")

        # 結果を図にして表示
        fig, ax = plt.subplots()
        im  = ax.imshow(levelset, vmin=np.min(levelset), vmax=np.max(levelset))
        cbar = fig.colorbar(im)
        cbar.set_label("Distance From phi = 0", fontsize=10)
        # plt.axis((self.min_x, self.max_x, self.min_y, self.max_y))
        levels = np.arange(0,0.01,1e-3)
        ctr = ax.contour(levelset, levels)#これを何回かごとに保存する．
        ax.clabel(ctr, levels, inline=1)
        plt.show()

        # 計算結果をcsvファイルに保存．（オプション）
        if settings["output_initallevelset"] and settings["mode"]=="geometry":
            output_filename = "init_levelset.csv" # + setting_filename.replace(".json",".csv")
            description  =f"settings : {setting_filename}\n" + "geometry : " + ", ".join([input_geometry["filename"] for input_geometry in settings["geometry"]["geometry"]])
            np.savetxt(output_filename, levelset, header=description, fmt='%.5f', delimiter=",")

        return levelset, A_p_init, l_p_init

    def culc_levelset_t_evo(self, levelset_old, rdot, del_t=0.001):
        self.r_arr = np.append(self.r_arr, rdot*del_t)
        levelset_new = levelset_old - rdot*del_t
        return levelset_new

    def culc_max_r(self, levelset_fin_filename):
        levelset_fin = np.loadtxt(levelset_fin_filename, delimiter=",", dtype=float, encoding='utf-8')
        delta_x = 0.1
        delta_y = 0.1 
        distance = distance_transform_edt(levelset_fin < 0, sampling=[delta_x, delta_y])
        max_r = max(distance)#ちがう
        return max_r

if __name__=='__main__':
    geom = FuelGeometry()
    levelset, A_p, l_p = geom.culc_initial_levelset("geometry_settings.jsonc")
    print(f"A_p:{A_p}\nl_p:{l_p}")
