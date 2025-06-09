import datetime
import numpy as np


class Satellite:
    def __init__(self, SatelliteName, Time, SatelliteClockCorrect, SatelliteObservation):
        # 定义基本的卫星信息
        self.SatelliteName = SatelliteName  # 卫星名称
        self.Time = Time  # 轨道参数参考时刻toe
        self.GpsSeconds = ctime2gps(Time)  # 将toe转换为gps周+周内秒形式
        self.SatelliteClockCorrect = SatelliteClockCorrect  # 卫星钟差参数
        self.SatelliteObservation = SatelliteObservation  # 卫星位置计算参数
        self.Delta_T = 0  # 卫星钟差改正数
        self.X = 0.0
        self.Y = 0.0
        self.Z = 0.0

    # 计算卫星观测时刻的位置函数
    def InitPositionOfSat(self, ObsTime):
        [sat_x, sat_y, sat_z, delta_t] = calculate_pos_of_sat(
            self.SatelliteObservation, self.SatelliteClockCorrect, ObsTime, self.Time
        )
        # 将得到的结果存入对象中
        self.X = sat_x
        self.Y = sat_y
        self.Z = sat_z
        self.Delta_T = delta_t

    # 获取卫星坐标的方法
    def GetPositionOfSat(self):
        return self.X, self.Y, self.Z


# 将时间转换为gps周+周内秒形式函数
def ctime2gps(time):
    gps_epoch = datetime.datetime(1980, 1, 6, 0, 0, 0)
    given_time = datetime.datetime(time[0], time[1], time[2], time[3], time[4], int(time[5]))
    time_diff = given_time - gps_epoch
    total_seconds = time_diff.total_seconds()
    gps_weeks, gps_seconds = divmod(total_seconds, 24 * 60 * 60 * 7)
    return gps_weeks, gps_seconds


# 计算卫星位置的静态方法
def calculate_pos_of_sat(matrix, ClockCorrect, ObsTime, RefTime):
    # 计算卫星钟差改正数
    gps_obs_weeks, gps_obs_sec = ctime2gps(ObsTime)
    gps_ref_weeks, gps_ref_sec = ctime2gps(RefTime)

    delta_t = (
            ClockCorrect[0]
            + ClockCorrect[1] * (gps_ref_sec - gps_obs_sec)
            + ClockCorrect[2] * pow(gps_ref_sec - gps_obs_sec, 2)
    )

    # 计算规化时刻
    t = gps_obs_sec - delta_t
    tk = t - matrix[2, 0]

    # 计算平近点角
    GM = 3.986005e14
    a = pow(matrix[1, 3], 2)  # 轨道长半轴
    n = np.sqrt(GM / pow(a, 3)) + matrix[0, 2]  # 平均角速度
    M = matrix[0, 3] + n * tk  # 平近点角
    M = M % (2 * np.pi)  # 将结果规范到[0, 2pi)，便于后续计算

    # 迭代求解偏近点角
    E_old = M
    for i in range(10):
        Ek = M + matrix[1, 1] * np.sin(E_old)
        if abs(Ek - E_old) < 1e-10:
            break
        else:
            E_old = Ek

    # 计算真近点角
    Vk = np.arctan(2 * np.sqrt(1 - pow(matrix[1, 1], 2)) * np.sin(Ek)
                   / (np.cos(Ek) - matrix[1, 1]))
    cos_vk = np.cos(Vk)
    sin_vk = np.sin(Vk)

    # 计算升交点角距
    phi_vk = Vk + matrix[3, 2]
    phi_vk = phi_vk % (2 * np.pi)   # 将结果规范到[0, 2pi)，便于后续计算

    # 计算摄动修正后的幅角、距离和轨道倾角
    delta_u = matrix[1, 0] * np.cos(2 * phi_vk) + matrix[1, 2] * np.sin(2 * phi_vk)    # 升交点角距改正项
    delta_r = matrix[3, 1] * np.cos(2 * phi_vk) + matrix[0, 1] * np.sin(2 * phi_vk)    # 轨道半径改正项
    delta_i = matrix[2, 1] * np.cos(2 * phi_vk) + matrix[2, 3] * np.sin(2 * phi_vk)    # 轨道倾角改正项
    uk = phi_vk + delta_u   # 升交点角距
    rk = a * (1 - matrix[1, 1] * np.cos(Ek)) + delta_r    # 轨道半径
    ik = matrix[3, 0] + delta_t + matrix[4, 0] * tk    # 轨道倾角

    # 计算轨道平面下的坐标
    xk = rk * np.cos(uk)
    yk = rk * np.sin(uk)

    # 地球自转角速度rad_v
    rad_v = 7.2921151467e-5
    # 升交点赤经
    omega_k = matrix[2, 2] + (matrix[3, 3] - rad_v) * tk - rad_v * matrix[2, 0]

    # 将轨道平面坐标旋转至地心地固坐标系，卫星的ECEF坐标
    Xk = xk * np.cos(omega_k) - yk * np.cos(ik) * np.sin(omega_k)
    Yk = xk * np.sin(omega_k) + yk * np.cos(ik) * np.cos(omega_k)
    Zk = yk * np.sin(ik)

    # 考虑相对论效应后的钟差
    F = -4.442807633e-10
    delta_t += F * matrix[1, 1] * np.sin(Ek) * matrix[1,3] - matrix[5, 2]

    return Xk, Yk, Zk, delta_t
