from satellite import Satellite
import numpy as np


class ReadFile:
    # 定义类变量

    # 粗略的接收机位置
    ApproxPos = []

    # O文件的所有行数据
    OLines = None

    # N文件的所有行数据
    NLines = None

    # O文件"END OF HEADER"所在行
    OHeaderLastLine = 0

    # N文件"END OF HEADER"所在行
    NHeaderLastLine = 0

    # 类的初始化函数
    def __init__(self, File):

        # 文件路径
        self.OFilePath = File[0]
        self.NFilePath = File[1]

        # 通过实例函数读取N、O文件的数据
        ReadFile.OLines = self.ReadOFile()
        ReadFile.NLines = self.ReadNFile()

        # 对O,N文件进行预处理
        ReadFile.OHeaderLastLine = ReadFile.PreprocessOFile(ReadFile.OLines)
        ReadFile.NHeaderLastLine = ReadFile.PreprocessNFile(ReadFile.NLines)

        # 存储卫星对象
        self.Satellites = []

        # 存储卫星名称、参考时间、卫星钟差改正、卫星轨道参数
        self.PosName = []
        self.Time = []
        self.SatelliteClockCorrect = []
        self.SatelliteObservation = []

    # 类方法的定义，用于获取各类变量的值
    @classmethod
    def GetApproxPos(cls):
        return cls.ApproxPos

    @classmethod
    def GetOLines(cls):
        return cls.OLines

    @classmethod
    def GetNLines(cls):
        return cls.NLines

    @classmethod
    def GetOHeaderLastLine(cls):
        return cls.OHeaderLastLine

    @classmethod
    def GetNHeaderLastLine(cls):
        return cls.NHeaderLastLine

    # 文件读取方法
    def ReadOFile(self):
        with open(self.OFilePath, 'r') as file:
            lines = file.readlines()
            ReadFile.OLines = lines
            return lines

    def ReadNFile(self):
        with open(self.NFilePath, 'r') as file:
            lines = file.readlines()
            ReadFile.NLines = lines
            return lines

    # O文件预处理函数，由于需要获取粗略坐标，所以需要传入O文件行数
    @staticmethod
    def PreprocessOFile(lines):
        ApproxPosComment = "APPROX POSITION XYZ"
        for i, line in enumerate(lines, start=1):
            if ApproxPosComment in line:
                approx_x = float(line[0:15].strip())
                approx_y = float(line[15:28].strip())
                approx_z = float(line[29:42].strip())
                ReadFile.ApproxPos.append(approx_x)
                ReadFile.ApproxPos.append(approx_y)
                ReadFile.ApproxPos.append(approx_z)
        ObsTargetString = "END OF HEADER"
        ObsHeaderLine = 0
        for i, line in enumerate(lines, start=1):
            if ObsTargetString in line:
                ObsHeaderLine = i
                break
        return ObsHeaderLine

    # N文件预处理函数
    @staticmethod
    def PreprocessNFile(lines):
        NavTargetString = "END OF HEADER"
        NavHeaderLine = 0
        for i, line in enumerate(lines, start=1):
            if NavTargetString in line:
                NavHeaderLine = i
                break
        return NavHeaderLine

    # 计算卫星参考时刻toe位置函数
    def CalculateSatellites(self):
        # 从N文件头部结束以后的第一个数据行开始，每8行为一个卫星的数据块，因此以8为步长循环处理
        for i in range(ReadFile.NHeaderLastLine, len(ReadFile.NLines) - 7, 8):
            # 每颗卫星第一行数据
            line = ReadFile.NLines[i]
            num = line[0:2].strip()
            # 存储时间年月日时分秒，这里的时间实际上是toc（钟差参数参考时刻）
            # 在现代大部分RINEX导航文件中toc和toe通常是一样的
            # 为保险起见将钟差参考时刻转换为周内秒形式后发现和toe（轨道参数参考时刻）相同
            # 因此之后卫星位置计算部分过程中直接使用了time数据表示toe
            time = [None] * 6
            time[0] = int((line[3:5]).strip()) + 2000
            time[1] = int((line[6:8]).strip())
            time[2] = int((line[9:11]).strip())
            time[3] = int((line[12:14]).strip())
            time[4] = int((line[15:17]).strip())
            time[5] = float((line[18:22]).strip())

            # 读取卫星钟差改正参数
            time_change = []
            # 卫星钟差af0
            a = float(line[22:37].strip()) * pow(10, int(line[38:41].strip()))
            time_change.append(a)
            # 卫星钟偏af1
            b = float(line[42:56].strip()) * pow(10, int(line[57:60]))
            time_change.append(b)
            # 卫星钟偏移af2
            c = float(line[60:75].strip()) * pow(10, int(line[76:79].strip()))
            time_change.append(c)

            self.SatelliteClockCorrect.append(time_change)

            # 读取卫星位置计算的参数，在一个卫星8行的数据块中，第一行 为和钟差有关的参数
            # 第八行为信息发射时间和拟合区间，在本次程序中不会用到，因此只读取其中6行4列的数据
            rows = 6
            cols = 4
            matrix = np.zeros((rows, cols))

            for j in range(0, rows):
                for k in range(0, cols):
                    matrix[j][k] = float(ReadFile.NLines[i + 1 + j][3 + 19 * k: 18 + 19 * k]
                                         ) * pow(10, int(ReadFile.NLines[i + 1 + j][19 + 19 * k:22 + 19 * k]))
            self.SatelliteObservation.append(matrix)

            # 保存卫星的名称
            SatelliteName = "G" + str(num)

            # 根据对应数据创建相应的卫星
            satellite = Satellite(SatelliteName, time, time_change, matrix)

            # 储存对应的卫星信息
            self.Time.append(time)
            self.PosName.append(SatelliteName)
            self.Satellites.append(satellite)
