from readfile import ReadFile
from satellite import Satellite
import datetime
import numpy as np


def CalculateTimeDifference(Time1, Time2):
    time1 = datetime.datetime(
        Time1[0], Time1[1], Time1[2], Time1[3], Time1[4], int(Time1[5])
    )
    time2 = datetime.datetime(
        Time2[0], Time2[1], Time2[2], Time2[3], Time2[4], int(Time2[5])
    )
    diff = time2 - time1
    seconds = diff.total_seconds()
    return seconds


def eecf2pos(ApproxPos):
    # WGS-84参考椭球的数据
    a = 637813.0  # 长半轴
    f = 1 / 298.257223563  # 扁率
    e2 = f * (2 - f)  # 偏心率平方

    Pos = [None] * 3

    # 求出平面距离
    L = np.sqrt(ApproxPos[0] * ApproxPos[0] + ApproxPos[1] * ApproxPos[1])
    # 初始纬度
    p0 = np.arctan(ApproxPos[2] / L)
    pi = 0.0
    # 迭代开始
    while 1:
        N = a / np.sqrt(1 - e2 * np.sin(p0) * np.sin(p0))
        pi = np.arctan((ApproxPos[2] + e2 * N * np.sin(p0)))
        if abs(pi - p0) < 1e-4:
            break
        else:
            p0 = pi

    # 计算高程h和经度lambda
    h = L / np.cos(pi) - N
    lam = np.arctan(ApproxPos[0] / ApproxPos[1])
    Pos.append(lam)
    Pos.append(pi)
    Pos.append(h)
    return Pos


def SolutionLeastSquares(ObsPseudorange, SatelliteXYZ, ApproxPos):
    # 首先需要判断一下匹配上的伪距个数是否大于4，由于进行最小二乘时有四个未知数，所以在求解时最少需要四个观测方程才能求解
    t = 0
    # 获取到的卫星数目
    SizeSatelliteXYZ = len(SatelliteXYZ)

    for i in range(SizeSatelliteXYZ):
        t = t + SatelliteXYZ[i][4]
    if t < 4:
        print("无法进行平差计算")
    else:
        dtr = 0.0
        # 进行平差迭代计算,i表示迭代次数
        for i in range(10):
            B = []
            L = []
            for k in range(SizeSatelliteXYZ):
                # 跳过没有匹配上的卫星
                if SatelliteXYZ[k][4] == 0:
                    continue

                # 结束标志
                if k == SizeSatelliteXYZ - 1:
                    break

                # 卫星坐标
                Xs, Ys, Zs = SatelliteXYZ[k][0:3]
                print(Xs, Ys, Zs)

                # 计算接收机到卫星的距离
                dx = Xs - ApproxPos[0]
                dy = Ys - ApproxPos[1]
                dz = Zs - ApproxPos[2]
                P0 = np.sqrt(dx * dx + dy * dy + dz * dz)

                # 计算Bj矩阵的每一行数据
                TemB3 = 1
                c = 299792458.0
                TemB0 = -dx / P0
                TemB1 = -dy / P0
                TemB2 = -dz / P0

                TemBRol = [TemB0, TemB1, TemB2, TemB3]
                B.append(TemBRol)

                # 计算L矩阵每一行的数据
                TemLRol = ObsPseudorange[k] - (P0 - c * SatelliteXYZ[k][3] + c * dtr)
                L.append(TemLRol)
                print(f"k={k}, 伪距={ObsPseudorange[k]}, P0={P0}, dt_sat={SatelliteXYZ[k][3]}, dtr={dtr}, L={TemLRol}")

            # 求解坐标改正量
            ArrayB = np.array(B)
            ArrayL = np.array(L).reshape(-1, 1)

            # 进行间接平差公式进行间接平差计算
            x, residuals, rank, s = np.linalg.lstsq(ArrayB, ArrayL, rcond=None)

            print(x)

            ApproxPos[0] = ApproxPos[0] + x[0, 0]
            ApproxPos[1] = ApproxPos[1] + x[1, 0]
            ApproxPos[2] = ApproxPos[2] + x[2, 0]
            dtr += x[3, 0]

            if np.linalg.norm(x[:3]) < 1e-4:
                break

        # 平差后的坐标存储
        PosXYZ = np.array(ApproxPos)
        print(
            "平差后的X坐标:",
            PosXYZ[0],
            "平差后的Y坐标:",
            PosXYZ[1],
            "平差后的Z坐标:",
            PosXYZ[2],
        )
    return 1


class Position:
    def __init__(self, SatelliteObservation, SatelliteName, Time, SatelliteClockCorrect):

        # N文件读取的所有卫星数据
        self.SatelliteObservation = SatelliteObservation
        self.SatelliteName = SatelliteName
        self.Time = Time
        self.SatelliteClockCorrect = SatelliteClockCorrect

        # 获取O文件的数据信息和头部结束行号
        self.Lines = ReadFile.OLines
        self.OHeaderLastLine = ReadFile.OHeaderLastLine

        # 获取观测值粗略坐标
        self.ApproxPos = ReadFile.ApproxPos

    # 观测数据和卫星数据进行匹配并计算测站位置
    def MatchObservationAndCalculate(self):
        # 在O文件中数据是以历元进行划分的，要进行ON文件卫星数据的匹配就需要知道每个历元数据行范围
        line = self.Lines[self.OHeaderLastLine]  # 历元数据块开始行数据
        NReadLine = self.OHeaderLastLine  # 历元数据块结束行行号

        # 进行循环匹配直到所有历元数据都被匹配完毕
        while line != "":
            # 读取观测时刻
            obs_time = [None] * 6
            obs_time[0] = int((line[1:3]).strip()) + 2000
            obs_time[1] = int((line[4:6]).strip())
            obs_time[2] = int((line[7:9]).strip())
            obs_time[3] = int((line[10:12]).strip())
            obs_time[4] = int((line[13:15]).strip())
            obs_time[5] = int((line[17:19]).strip())

            # 读取本历元观测到的卫星数
            num_sat = int(line[30:32])

            # 读取本历元中观测到的卫星名称
            # 由于O文件中一行最多只能存储12个卫星名称，因此判断一下卫星名称有几行
            if num_sat % 12 == 0:
                n = int(num_sat / 12)
            else:
                n = int(num_sat / 12) + 1

            PRN_str = ""
            for j in range(n):
                if j == 0:
                    PRN_str = PRN_str + self.Lines[NReadLine + j][32:69]
                else:
                    PRN_str = PRN_str + self.Lines[NReadLine + j][0:].strip()
            obs_sat_PRN = []
            for k in range(num_sat):
                obs_sat_PRN.append(PRN_str[k * 3: k * 3 + 3])

            NReadLine = NReadLine + n

            # 读取卫星伪距观测值
            ObsPseudorange = []

            for j in range(NReadLine, NReadLine + 2 * num_sat, 2):
                line = self.Lines[j]
                Pseudorange = [None] * 4
                Pseudorange[0] = line[34:46].strip()
                Pseudorange[1] = line[50:62].strip()
                Pseudorange[2] = line[54:66].strip()
                Pseudorange[3] = line[70:].strip()
                for i in range(4):
                    if Pseudorange[i] == "":
                        Pseudorange[i] = 0
                    else:
                        Pseudorange[i] = float(Pseudorange[i])
                ObsPseudorange.append(Pseudorange)

            # 直接计算当前位置下的接收机位置
            SatelliteXYZ, M_ObsPseudorange = self.MatchToSatellite(
                obs_time, obs_sat_PRN, ObsPseudorange
            )

            # 运用非线性最小二乘，平差计算地面坐标
            a = SolutionLeastSquares(
                M_ObsPseudorange, SatelliteXYZ, ReadFile.ApproxPos
            )

            # 更新读取的行数，和行的内容
            NReadLine = NReadLine + 2 * num_sat
            line = self.Lines[NReadLine]

    def MatchToSatellite(self, obs_time, obs_sat_PRN, ObsPseudorange):
        SatelliteXYZ = []
        M_ObsPseudorange = []
        # 遍历观测文件指定历元的卫星编号
        for index, SatPRN in enumerate(obs_sat_PRN):
            # 将观测文件中的卫星编号和导航文件中的卫星编号进行匹配

            # 存放五个元素，前四个分别为卫星坐标和卫星钟差，最后一个是是否匹配成功的标记
            TemXYZ = [None] * 5

            TimeDiff = []
            for index1, SatPRN1 in enumerate(self.SatelliteName):
                if SatPRN == SatPRN1:
                    # 计算两个时间差
                    TimeDiff.append(
                        CalculateTimeDifference(obs_time, self.Time[index1])
                    )
                else:
                    # 这里设置三天所有的秒数，保证足够大，找最小值的之后不找到就可以了
                    TimeDiff.append(2592000)

            # 判断是否匹配成功
            NotMatch = all(x == 2592000 for x in TimeDiff)
            if NotMatch:
                TemXYZ[0] = 0
                TemXYZ[1] = 0
                TemXYZ[2] = 0
                TemXYZ[3] = 0
                TemXYZ[4] = 0
                M_ObsPseudorange.append(0.0)
            else:
                # 寻找最小时间索引
                MinTime = min(TimeDiff)
                MinTimeIndex = TimeDiff.index(MinTime)
                print(self.SatelliteName[MinTimeIndex])
                # 获得这个卫星的信息
                satellite = Satellite(
                    self.SatelliteName[MinTimeIndex],
                    self.Time[MinTimeIndex],
                    self.SatelliteClockCorrect[MinTimeIndex],
                    self.SatelliteObservation[MinTimeIndex],
                )
                # 传入观测时间，根据初步改正的伪距观测值计算信号发射时刻从而计算卫星位置
                c = 299792458.0
                f_L1 = 1575.42 * pow(10, 6)
                f_L2 = 1227.60 * pow(10, 6)
                Lambda_L1 = c / f_L1
                Lambda_L2 = c / f_L2
                Gamma = pow((Lambda_L2 / Lambda_L1), 2)
                P1 = ObsPseudorange[index][0]
                P2 = ObsPseudorange[index][1]
                C1 = ObsPseudorange[index][2]
                C2 = ObsPseudorange[index][3]
                # 进行双频或者时单频改正
                if not P1 == 0:
                    P1_obs = P1
                elif not C1 == 0:
                    P1_obs = C1
                else:
                    TemXYZ = [0, 0, 0, 0, 0]
                    SatelliteXYZ.append(TemXYZ)
                    M_ObsPseudorange.append(0.0)
                    continue  # 缺失L1观测直接放弃本星
                if not P2 == 0:
                    P2_obs = P2
                elif not C2 == 0:
                    P2_obs = C2
                else:
                    P2_obs = 0  # 是否缺失L2观测作为判断进行双频还是单频改正的依据
                if not P2_obs == 0:
                    Pc = (Gamma * P1_obs - P2_obs) / (Gamma - 1.0)  # 双频改正
                else:
                    P1_P2 = (1.0 - Gamma) * satellite.SatelliteObservation[2, 5]  # 单频改正
                    Pc = P1_obs - P1_P2 * (1.0 - Gamma)
                M_ObsPseudorange.append(Pc)  # 将初步改正后的伪距存入

                Obs_time = datetime.datetime(obs_time[0], obs_time[1], obs_time[2], obs_time[3], obs_time[4],
                                             int(obs_time[5]))
                Pro_time = Pc / c
                time_delta = datetime.timedelta(seconds=Pro_time)
                transmit_time = Obs_time - time_delta  # # 通过伪距计算传播时间，得到真正的信号发射时刻
                M_obs_time = [transmit_time.year, transmit_time.month, transmit_time.day,
                              transmit_time.hour, transmit_time.minute, transmit_time.second
                              ]

                satellite.InitPositionOfSat(M_obs_time)
                TemXYZ[0] = satellite.X
                TemXYZ[1] = satellite.Y
                TemXYZ[2] = satellite.Z
                TemXYZ[3] = satellite.Delta_T
                TemXYZ[4] = 1
            SatelliteXYZ.append(TemXYZ)
        return SatelliteXYZ, M_ObsPseudorange
