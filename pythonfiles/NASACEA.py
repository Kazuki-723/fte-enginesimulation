import subprocess
import os
import re

def CEA(Pc, OF,epsilon):
    path = 'data\\yoshidaexcel.inp'

    #単位変換、小数点以下を丸める
    #P_cはbarrで入れる
    Pc = 10 * Pc
    Pc = round(Pc,3)
    OF = round(OF,3)
    #P_inj = round(P_c/1.013,3)



    with open(path, mode = "w") as f:
        f.write("problem    o/f= ")
        f.write(str(OF))
        f.write("\n")

        f.write("      rocket sup,ae/at = ")
        f.write(str(epsilon))
        #f.write(" nfz = 2")
        f.write("\n")

        f.write("  p,bar=")
        f.write(str(Pc))
        f.write("\n")

        #todo oxのデータとfuelのデータを外部入力可能にする
        f.write("react\n")
        f.write(" oxid=N2O wt=100  t,k=290\n")
        f.write(" fuel=PMMA wt=100  t,k=290\n")
        f.write("    h,kj/mol=-468.3  C 5 O 2 H 8\n")

        #ABSのデータ
        # f.write(" fuel=acrylonitrile wt=43  t,k=290\n")
        # f.write("    h,kj/mol=172.6  C 3 N 1 H 3\n")
        # f.write(" fuel=butadiene wt=47  t,k=290\n")
        # f.write("    h,kj/mol=108.8  C 3 H 6\n")
        # f.write(" fuel=styrene wt=10  t,k=290\n")
        # f.write("    h,kj/mol=146.9  C 8 H 8\n")

        f.write("end")

    #CEAを走らせる
    os.chdir('data')
    with subprocess.Popen('FCEA2m.exe', stdin=subprocess.PIPE, encoding="ASCII", stdout = subprocess.DEVNULL) as proc:
        proc.stdin.write("yoshidaexcel\n")
    os.chdir('..')

    #output抽出
    path = 'data\\yoshidaexcel.out'
    outputlist = []

    with open(path) as f:
        #各行を収納した配列の生成
        for s_line in f:
            outputlist.append(repr(s_line))
        
    #改行コードのみの行を削除
    outputlist_new = [a for a in outputlist if a != "'\\n'"]

    #必要なパラメータを取り出す C_F,Cstar,gamma,T_c
    #Ae/Atがepsilonを示してる

    for i in range(len(outputlist_new)):
        if 'GAMMA' in outputlist_new[i]:
            gamma_num = re.findall(r"\d*\.\d+",outputlist_new[i])
        elif 'CSTAR' in outputlist_new[i]:
            Cstar_num = re.findall(r"\d*\.\d+",outputlist_new[i])
        elif 'CF' in outputlist_new[i]:
            CF_num = re.findall(r"\d*\.\d+",outputlist_new[i])
        elif 'T, K' in outputlist_new[i]:
            T_c_num = re.findall(r"\d*\.\d+",outputlist_new[i])
        elif 'M, (1/n)' in outputlist_new[i]:
            Mole_num = re.findall(r"\d*\.\d+",outputlist_new[i])
        # elif 'Ae/At' in outputlist_new[i]:
        #     epsilon_num = re.findall(r"\d*\.\d+",outputlist_new[i])
        elif 'P, BAR' in outputlist_new[i]:
            Pt_num = re.findall(r"\d*\.\d+",outputlist_new[i])
        elif 'MACH NUMBER' in outputlist_new[i]:
            Mach_num = re.findall(r"\d*\.\d+",outputlist_new[i])
    
    #0:chamber 1:throat 2:exit
    #rocket parametorsはchamberの値がないので、0:throat 1:exit 
    gamma_new = float(gamma_num[2])#出口比熱比
    Cstar_new = float(Cstar_num[1])#出口Cstar
    CF_new = float(CF_num[1] )#出口CF
    T_c_new = float(T_c_num[0]) #燃焼室温度
    T_t_new = float(T_c_num[1]) #スロート温度
    T_e_new = float(T_c_num[2]) #出口温度
    Mole_new = float(Mole_num[2])
    # epsilon_new = float(epsilon_num[1]) #出口開口比
    Pthroat_new = float(Pt_num[1])/10 #MPaに直す,スロート圧
    Mach_new = float(Mach_num[2])#出口速度マッハ数
    Pe_new = float(Pt_num[2])/10#出口圧

    return gamma_new, Cstar_new, CF_new, T_c_new, T_t_new, T_e_new, Mole_new, Pthroat_new, Pe_new, Mach_new