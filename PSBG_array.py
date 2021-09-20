###########################################################################
####  run (1) SiEPIC_EBeam_Library.lym  (2) MyLib.lym  (3) this file ############################
###########################################################################
import pya
from SiEPIC.scripts import connect_cell,connect_pins_with_waveguide,path_to_waveguide
from SiEPIC.utils import get_layout_variables
import numpy as np

def flatten(t):
    return [item for sublist in t for item in sublist]

class DrawLayout(object):
    def __init__(self):
        # prepare
        self.TECHNOLOGY, self.lv, self.ly, self.cell = get_layout_variables()
        self.TECHNOLOGY = get_technology_by_name('EBeam')
        # LayerInfo
        self.LayerSi = pya.LayerInfo(1, 0)
        self.LayerSiN = self.cell.layout().layer(self.LayerSi)
        self.TextLayerN = self.cell.layout().layer(pya.LayerInfo(10, 0))
        self.FloorPlan = self.cell.layout().layer(pya.LayerInfo(99, 0))

        self.tc = self.ly.top_cells()[0] # top cell
        self.ly.prune_subcells(self.tc.cell_index(), 10)


    def FP(self, x_span, y_span, x_c, y_c):
        # top cell

        box = pya.Box(x_c-x_span/2, y_c-y_span/2, x_c + x_span/2, y_c + y_span/2)
        # box = pya.Box(0, 0, 1000, 2000)
        self.tc.shapes(self.FloorPlan).insert(box)
        return self.tc

    def SWG_PSBG_cell(self, psbg_name, dir_path, filenames, param_psbg1, cell_x, cell_y):

        cell_psbg = self.ly.create_cell(psbg_name)
        cell_psbg_core = self.ly.create_cell("psbg_core")
        # insert psbg/instance

        param_psbg1["Path_1"] = dir_path + filenames[0]
        param_psbg1["Path_2"] = dir_path + filenames[1]

        psbg1 = self.ly.create_cell("Chirped_SWG_BG_wing", "MyLib", param_psbg1)
        t_psbg1 = pya.Trans.from_s('r0 0,0')
        inst_psbg1 = cell_psbg_core.insert(pya.DCellInstArray(psbg1.cell_index(), t_psbg1))

        #insert taper
        param_taper1 = {"wavelength": 1550,
                "fishbone": False,
                "length": 25,
                "taper_fraction": 0.8,
                "period_strip": param_psbg1["SWG_period"],
                "period_swg": param_psbg1["SWG_period"],
                "wg_width_strip": param_psbg1["wg_width"],
                "wg_width_swg": param_psbg1["wg_width"],
                "wg_width_taper": 0.06,
                "duty_strip": param_psbg1["leng_swg"]/param_psbg1["SWG_period"],
                "duty_swg": param_psbg1["leng_swg"]/param_psbg1["SWG_period"]}
        taper1 = self.ly.create_cell("Waveguide_SWG_to_Strip", "EBeam_Beta", param_taper1)
        # insert wg
        param_wg1 = {"wg_width1": 0.5,
              "wg_width2": 0.4,
              "wg_length": 10}
        wg1 = self.ly.create_cell("ebeam_taper_te1550", "EBeam", param_wg1)
        # insert Y branch
        yb1 = self.ly.create_cell("ebeam_y_1550", "EBeam")

        # connect components from the left side of PSBG
        inst_tp = connect_cell(inst_psbg1, 'pin1', taper1, 'pin1')
        inst_wtp = connect_cell(inst_tp, 'pin2', wg1, 'pin2')
        inst_ywtp = connect_cell(inst_wtp, 'pin1', yb1, 'opt1')
        # connect components from the right side of PSBG
        inst_pt = connect_cell(inst_psbg1, "pin2", taper1, "pin1")
        inst_ptw = connect_cell(inst_pt, "pin2", wg1, "pin2")  # the pin label of inst_ptw (instance) is same to wg1


        #insert vgc
        y_vgc1 = 60
        vgc_space = 127
        self.ly.read("C:\\Users\\ZXY\\KLayout\\tech\\NanoSOI_PDK_v5\\tech\\libraries\\ANT_PDK_Silicon_v5.1.gds")
        vgc1 =self.ly.cell("GratingCoupler_TE_Oxide_8degrees")

        t_vgc1 = pya.Trans.from_s('r90 -100,'+str(y_vgc1))
        inst_vgc1 = cell_psbg.insert(pya.DCellInstArray(vgc1.cell_index(),  t_vgc1))
        t_vgc2 = pya.Trans.from_s('r90 -100,'+str(y_vgc1-vgc_space))
        inst_vgc2 = cell_psbg.insert(pya.DCellInstArray(vgc1.cell_index(),  t_vgc2))
        t_vgc3 = pya.Trans.from_s('r90 -100,'+str(y_vgc1-2*vgc_space))
        inst_vgc3 = cell_psbg.insert(pya.DCellInstArray(vgc1.cell_index(),  t_vgc3))

        # Waveguides:
        inst_vgc1_pin1 = inst_vgc1.pinPoint('pin1')
        inst_vgc2_pin1 = inst_vgc2.pinPoint('pin1')
        inst_vgc3_pin1 = inst_vgc3.pinPoint('pin1')
        inst_ptw_pin1 = inst_ptw.pinPoint('pin1')
        inst_ywtp_opt3 = inst_ywtp.pinPoint('opt3')
        inst_ywtp_opt2 = inst_ywtp.pinPoint('opt2')

        # add path
        wg_width = 0.5
        path1 = Path([inst_vgc1_pin1, Point(inst_vgc1_pin1.x+1000*20, inst_vgc1_pin1.y),Point(inst_vgc1_pin1.x+1000*20, inst_ywtp_opt3.y),inst_ywtp_opt3], wg_width)
        cell_psbg.shapes(self.LayerSiN).insert(path1)

        path2 = Path([inst_vgc2_pin1, Point(inst_vgc2_pin1.x+1000*20, inst_vgc2_pin1.y),Point(inst_vgc1_pin1.x+1000*20, inst_ywtp_opt2.y),inst_ywtp_opt2], wg_width)
        cell_psbg.shapes(self.LayerSiN).insert(path2)

        path3 = Path([inst_vgc3_pin1, Point(inst_ptw_pin1.x+1000*20, inst_vgc3_pin1.y),Point(inst_ptw_pin1.x+1000*20, inst_ptw_pin1.y),inst_ptw_pin1], wg_width)
        cell_psbg.shapes(self.LayerSiN).insert(path3)
        # path to wg
        wg_width_um = 0.5
        param_p2wg = {'width': wg_width_um, 'adiabatic': True, 'radius': 20.0, 'bezier': 0.2, 'wgs': [{'width': 0.5, 'layer': 'Si', 'offset': 0.0}]}
        path_to_waveguide(params = param_p2wg, cell=cell_psbg, snap=True, verbose=False)
        # add Label
        dbu = self.ly.dbu
        t_text = pya.Trans(pya.Trans.R0, (-150)/dbu, -150/dbu)
        text = pya.Text (psbg_name, t_text)
        shape = cell_psbg.shapes(self.TextLayerN).insert(text)
        shape.text_size = 25/dbu

        t_psbg_core = pya.Trans.from_s('r0 0,0')
        cell_psbg.insert(pya.DCellInstArray(cell_psbg_core.cell_index(), t_psbg_core))
        t_cell_psbg = pya.Trans.from_s('r0 '+str(cell_x)+','+str(cell_y))
        self.tc.insert(pya.DCellInstArray(cell_psbg.cell_index(),  t_cell_psbg))

        return self.tc

    def SWG_PSBG_array(self, psbg_name_list, dir_path, filename_list, param_psbg1, cell_x, cell_y, n_x, n_y, xspace, yspace):
        i = 0
        number = len(psbg_name_list)
        if n_x*n_y != number:
            print('The number of PSBGs (n_x*n_y) doesnot match the psbg_name_list')
        else:
            for nx in range(n_x):
                x = cell_x + nx*xspace
                for ny in range(n_y):
                    y = cell_y - ny*yspace
                    filenames = filename_list[i]
                    psbg_name = psbg_name_list[i]
                    self.SWG_PSBG_cell(psbg_name, dir_path, filenames, param_psbg1, x, y)
                    i = i+1

        return self.tc

    def Tilted_SWGBG_cell(self, tswgbg_name, param_tswgbg1, cell_x, cell_y, vgc_mode):

        cell_tswgbg1 = self.ly.create_cell(tswgbg_name)
        cell_tswgbg_core = self.ly.create_cell("tswgbg_core")
        # insert psbg/instance

        tswgbg1 = self.ly.create_cell("Tilted_SWG_BG_wing", "MyLib", param_tswgbg1)
        t_tswgbg1 = pya.Trans.from_s('r0 0,0')
        inst_tswgbg1 = cell_tswgbg_core.insert(pya.DCellInstArray(tswgbg1.cell_index(), t_tswgbg1))


        #insert taper
        param_taper1 = {"N_period": 25,
                "SWG_period": param_tswgbg1["SWG_period"],
                "wg_width":param_tswgbg1["wg_width"],
                "duty": param_tswgbg1["duty"],
                "theta": param_tswgbg1["theta"],
                "period_swg": param_tswgbg1["SWG_period"],
                "wg_width_strip": param_tswgbg1["wg_width"],
                "wg_width_swg": param_tswgbg1["wg_width"],
                "wg_width_taper": 0.06,
                "perc": 0.8}

        taper1 = self.ly.create_cell("Tilted_SWG_to_strip_waveguide", "MyLib", param_taper1)
        # insert wg
        param_wg1 = {"wg_width1": 0.5,
              "wg_width2": param_tswgbg1["wg_width"],
              "wg_length": 10}
        wg1 = self.ly.create_cell("ebeam_taper_te1550", "EBeam", param_wg1)
        # insert Y branch
        yb1 = self.ly.create_cell("ebeam_y_1550", "EBeam")

        # connect components from the left side of PSBG
        inst_tp = connect_cell(inst_tswgbg1, 'pin1', taper1, 'pin2')
        inst_wtp = connect_cell(inst_tp, 'pin1', wg1, 'pin2')
        inst_ywtp = connect_cell(inst_wtp, 'pin1', yb1, 'opt1')
        # connect components from the right side of PSBG
        inst_pt = connect_cell(inst_tswgbg1, "pin2", taper1, "pin2")
        inst_ptw = connect_cell(inst_pt, "pin1", wg1, "pin2")  # the pin label of inst_ptw (instance) is same to wg1


        #insert vgc

        y_vgc1 = 60
        vgc_space = 127
        self.ly.read("C:\\Users\\ZXY\\KLayout\\tech\\NanoSOI_PDK_v5\\tech\\libraries\\ANT_PDK_Silicon_v5.1.gds")
        if vgc_mode == 'TE':
            vgc1 =self.ly.cell("GratingCoupler_TE_Oxide_8degrees")

        elif vgc_mode == 'TM':
            vgc1 =self.ly.cell("GratingCoupler_TM_Oxide_8degrees")
        else:
            print("Please enter 'TE' or 'TM' at vgc_mode")


        t_vgc1 = pya.Trans.from_s('r90 -100,'+str(y_vgc1))
        inst_vgc1 = cell_tswgbg1.insert(pya.DCellInstArray(vgc1.cell_index(),  t_vgc1))
        t_vgc2 = pya.Trans.from_s('r90 -100,'+str(y_vgc1-vgc_space))
        inst_vgc2 = cell_tswgbg1.insert(pya.DCellInstArray(vgc1.cell_index(),  t_vgc2))
        t_vgc3 = pya.Trans.from_s('r90 -100,'+str(y_vgc1-2*vgc_space))
        inst_vgc3 = cell_tswgbg1.insert(pya.DCellInstArray(vgc1.cell_index(),  t_vgc3))

        # Waveguides:
        inst_vgc1_pin1 = inst_vgc1.pinPoint('pin1')
        inst_vgc2_pin1 = inst_vgc2.pinPoint('pin1')
        inst_vgc3_pin1 = inst_vgc3.pinPoint('pin1')
        inst_ptw_pin1 = inst_ptw.pinPoint('pin1')
        inst_ywtp_opt3 = inst_ywtp.pinPoint('opt3')
        inst_ywtp_opt2 = inst_ywtp.pinPoint('opt2')

        # add path
        wg_width = 0.5
        path1 = Path([inst_vgc1_pin1, Point(inst_vgc1_pin1.x+1000*20, inst_vgc1_pin1.y),Point(inst_vgc1_pin1.x+1000*20, inst_ywtp_opt3.y),inst_ywtp_opt3], wg_width)
        cell_tswgbg1.shapes(self.LayerSiN).insert(path1)

        path2 = Path([inst_vgc2_pin1, Point(inst_vgc2_pin1.x+1000*20, inst_vgc2_pin1.y),Point(inst_vgc1_pin1.x+1000*20, inst_ywtp_opt2.y),inst_ywtp_opt2], wg_width)
        cell_tswgbg1.shapes(self.LayerSiN).insert(path2)

        path3 = Path([inst_vgc3_pin1, Point(inst_ptw_pin1.x+1000*20, inst_vgc3_pin1.y),Point(inst_ptw_pin1.x+1000*20, inst_ptw_pin1.y),inst_ptw_pin1], wg_width)
        cell_tswgbg1.shapes(self.LayerSiN).insert(path3)
        # path to wg
        wg_width_um = 0.5
        param_p2wg = {'width': wg_width_um, 'adiabatic': True, 'radius': 20.0, 'bezier': 0.2, 'wgs': [{'width': 0.5, 'layer': 'Si', 'offset': 0.0}]}
        path_to_waveguide(params = param_p2wg, cell=cell_tswgbg1, snap=True, verbose=False)
        # add Label
        dbu = self.ly.dbu
        if vgc_mode =='TE':
            t_text = pya.Trans(pya.Trans.R0, (-150)/dbu, -150/dbu)
        else:
            t_text = pya.Trans(pya.Trans.M90, (280)/dbu, -150/dbu)
        text = pya.Text (tswgbg_name, t_text)
        shape = cell_tswgbg1.shapes(self.TextLayerN).insert(text)
        shape.text_size = 10/dbu

        t_tswgbg1_core = pya.Trans.from_s('r0 0,0')
        cell_tswgbg1.insert(pya.DCellInstArray(cell_tswgbg_core.cell_index(), t_tswgbg1_core))
        if vgc_mode == 'TE':
            t_cell_tswgbg1 = pya.Trans.from_s('r0 '+str(cell_x)+','+str(cell_y))
        else:
            t_cell_tswgbg1 = pya.Trans.new(pya.Trans.M90,(cell_x + 160), cell_y)#('r0,bool mirr = true,'+str(cell_x)+','+str(cell_y))
        self.tc.insert(pya.DCellInstArray(cell_tswgbg1.cell_index(),  t_cell_tswgbg1))

        return self.tc

    def Tilted_SWGBG_TETM_array(self, param_tswgbg, cell_x, cell_y, n_x, n_y, xspace, yspace):
        i = 0
        n_period = param_tswgbg["N_period"]
        swg_period = param_tswgbg["SWG_period"]
        wg_width = param_tswgbg["wg_width"]
        duty = param_tswgbg["duty"]
        theta = param_tswgbg["theta"]
        wing1_width = param_tswgbg["Wing1_width"]
        wing2_width = param_tswgbg["Wing2_width"]
        distance_2 = param_tswgbg["distance_2"]

        number = np.size(n_period)*np.size(swg_period)*np.size(wg_width)*np.size(duty)*np.size(theta)*np.size(wing1_width)*np.size(wing2_width)*np.size(distance_2)
        # if n_x*n_y != number:
        #     print('The number of tilted SWG BGs (n_x*n_y) doesnot match the number of parameters swept')
        # else:
        i_nx = 1
        i_ny = 1
        cell_x_i = cell_x
        cell_y_i = cell_y
        for i_N in range(np.size(n_period)):
            for i_S in range(np.size(swg_period)):
                for i_wg in range(np.size(wg_width)):
                    for i_duty in range(np.size(duty)):
                        for i_theta in range(np.size(theta)):
                            for i_W1 in range(np.size(wing1_width)):
                                for i_W2 in range(np.size(wing2_width)):
                                    for i_distance in range(np.size(distance_2)):
                                        tswgbg_name = 'TSWGBG_P'+str(swg_period[i_S]*1000)+'_W'+str(wg_width[i_wg])+'_DC'+str(duty[i_duty])+'_theta'+str(theta[i_theta])+'_Wing1W'+str(wing1_width[i_W1])+'_Wing2W'+str(wing2_width[i_W1])+'_gap'+str(distance_2[i_distance])
                                        param_tswgbg_i = {"N_period": n_period[i_N],
                                                         "SWG_period": swg_period[i_S],
                                                         "wg_width": wg_width[i_wg],
                                                         "duty": duty[i_duty],
                                                         "theta": theta[i_theta],
                                                         "Wing1_width": wing1_width[i_W1],
                                                         "Wing2_width": wing2_width[i_W2],
                                                         "distance_2": distance_2[i_distance]}
                                        self.Tilted_SWGBG_cell(tswgbg_name+'TE', param_tswgbg_i, cell_x_i, cell_y_i, 'TE')
                                        self.Tilted_SWGBG_cell(tswgbg_name+'TM', param_tswgbg_i, cell_x_i, cell_y_i-yspace, 'TM')
                                        i_ny = i_ny + 1
                                        if i_ny <= n_y:
                                            cell_y_i = cell_y_i - 2*yspace
                                        else:
                                            cell_x_i = cell_x_i + xspace
                                            cell_y_i = cell_y
                                            i_ny = 1

        return self.tc


####### prepare #########
# chip footprint
x_span = 9000.0 *1000#nm
y_span = 9000.0 *1000# nm
x_c = 0.0
y_c = 0.0

# first psbg cell
cell_x = round(-x_span/1000/2 + 150) #int
cell_y = round(y_span/1000/2 - 80)

# psbg array
n_x = 8# number of column
n_y = 1
xspace = 500 # micron  space in x direction
yspace = 300

group_space_y = 300

# psbg parameters

cell_x_P240 = cell_x #int
cell_y_P240 = cell_y

param_psbg_P240 =  {"N_period": 1000,
                "SWG_period": 0.240, # micron
                "wg_width": 0.4,
                "Wing1_width": 0.120,
                "Wing2_width": 0.120,
                "leng_swg": 0.120}

filename_list_P240 = [['PSBG_P240_sin1_gap_n1_c1500_BW3T.txt', 'PSBG_P240_sin1_gap_n2_c1500_BW3T.txt'],
             ['PSBG_P240_sin1_gap_n1_c1510_BW3T.txt', 'PSBG_P240_sin1_gap_n2_c1510_BW3T.txt'],
             ['PSBG_P240_sin1_gap_n1_c1520_BW3T.txt', 'PSBG_P240_sin1_gap_n2_c1520_BW3T.txt'],
             ['PSBG_P240_sin1_gap_n1_c1530_BW3T.txt', 'PSBG_P240_sin1_gap_n2_c1530_BW3T.txt'],
             ['PSBG_P240_sin1_gap_n1_c1540_BW3T.txt', 'PSBG_P240_sin1_gap_n2_c1540_BW3T.txt'],
             ['PSBG_P240_sin1_gap_n1_c1540_BW6T.txt', 'PSBG_P240_sin1_gap_n2_c1540_BW6T.txt'],
             ['PSBG_P240_sin1_gap_n1_c1530_BW6T.txt', 'PSBG_P240_sin1_gap_n2_c1530_BW6T.txt'],
             ['PSBG_P240_sin1_gap_n1_c1520_BW6T.txt', 'PSBG_P240_sin1_gap_n2_c1520_BW6T.txt']]

psbg_name_list_P240 = ['PSBG_c1500_BW3T_P240',
                    'PSBG_c1510_BW3T_P240',
                    'PSBG_c1520_BW3T_P240',
                    'PSBG_c1530_BW3T_P240',
                    'PSBG_c1540_BW3T_P240',
                    'PSBG_c1540_BW6T_P240',
                    'PSBG_c1530_BW6T_P240',
                    'PSBG_c1520_BW6T_P240']

# psbg parameters
cell_x_P245 = cell_x #int
cell_y_P245 = cell_y - group_space_y

param_psbg1 =  {"N_period": 1000,
                "SWG_period": 0.245, # micron
                "wg_width": 0.4,
                "Wing1_width": 0.122,
                "Wing2_width": 0.122,
                "leng_swg": 0.122}

# loading gap file infor
dir_path = 'F:\\SWG\\Code\\PSBG\\Layout_list\\gap_list_Apr6_2021_refab\\'
filename_list = [['PSBG_sin1_gap_n1_c1520_BW3T.txt', 'PSBG_sin1_gap_n2_c1520_BW3T.txt'],
             ['PSBG_sin1_gap_n1_c1530_BW3T.txt', 'PSBG_sin1_gap_n2_c1530_BW3T.txt'],
             ['PSBG_sin1_gap_n1_c1540_BW3T.txt', 'PSBG_sin1_gap_n2_c1540_BW3T.txt'],
             ['PSBG_sin1_gap_n1_c1550_BW3T.txt', 'PSBG_sin1_gap_n2_c1550_BW3T.txt'],
             ['PSBG_sin1_gap_n1_c1560_BW3T.txt', 'PSBG_sin1_gap_n2_c1560_BW3T.txt'],
             ['PSBG_P245_sin1_gap_n1_c1540_BW6T.txt', 'PSBG_P245_sin1_gap_n2_c1540_BW6T.txt'],
             ['PSBG_P245_sin1_gap_n1_c1550_BW6T.txt', 'PSBG_P245_sin1_gap_n2_c1550_BW6T.txt'],
             ['PSBG_P245_sin1_gap_n1_c1560_BW6T.txt', 'PSBG_P245_sin1_gap_n2_c1560_BW6T.txt']]

psbg_name_list = ['PSBG_c1520_BW3T_P245',
             'PSBG_c1530_BW3T_P245',
             'PSBG_c1540_BW3T_P245',
             'PSBG_c1550_BW3T_P245',
             'PSBG_c1560_BW3T_P245',
             'PSBG_c1540_BW6T_P245',
             'PSBG_c1550_BW6T_P245',
             'PSBG_c1560_BW6T_P245']

# psbg parameters
# first psbg cell
cell_x_P250 = cell_x #int
cell_y_P250 = cell_y - 2*group_space_y

param_psbg_P250 =  {"N_period": 1000,
                "SWG_period": 0.250, # micron
                "wg_width": 0.4,
                "Wing1_width": 0.125,
                "Wing2_width": 0.125,
                "leng_swg": 0.125}

filename_list_P250 = [['PSBG_P250_sin1_gap_n1_c1550_BW3T.txt', 'PSBG_P250_sin1_gap_n2_c1550_BW3T.txt'],
                     ['PSBG_P250_sin1_gap_n1_c1560_BW3T.txt', 'PSBG_P250_sin1_gap_n2_c1560_BW3T.txt'],
                     ['PSBG_P250_sin1_gap_n1_c1570_BW3T.txt', 'PSBG_P250_sin1_gap_n2_c1570_BW3T.txt'],
                     ['PSBG_P250_sin1_gap_n1_c1580_BW3T.txt', 'PSBG_P250_sin1_gap_n2_c1580_BW3T.txt'],
                     ['PSBG_P250_sin1_gap_n1_c1590_BW3T.txt', 'PSBG_P250_sin1_gap_n2_c1590_BW3T.txt'],
                     ['PSBG_P250_sin1_gap_n1_c1590_BW6T.txt', 'PSBG_P250_sin1_gap_n2_c1590_BW3T.txt'],
                     ['PSBG_P250_sin1_gap_n1_c1580_BW6T.txt', 'PSBG_P250_sin1_gap_n2_c1580_BW3T.txt'],
                     ['PSBG_P250_sin1_gap_n1_c1570_BW6T.txt', 'PSBG_P250_sin1_gap_n2_c1570_BW3T.txt']]

psbg_name_list_P250 = ['PSBG_c1550_BW3T_P250',
             'PSBG_c1560_BW3T_P250',
             'PSBG_c1570_BW3T_P250',
             'PSBG_c1580_BW3T_P250',
             'PSBG_c1590_BW3T_P250',
             'PSBG_c1590_BW6T_P250',
             'PSBG_c1580_BW6T_P250',
             'PSBG_c1570_BW6T_P250']

# psbg parameters
# first psbg cell
cell_x_P254 = cell_x + 9*xspace#int
cell_y_P254 = cell_y

param_psbg_P254 =  {"N_period": 1000,
                "SWG_period": 0.254, # micron
                "wg_width": 0.4,
                "Wing1_width": 0.127,
                "Wing2_width": 0.127,
                "leng_swg": 0.127}

filename_list_P254 = [['PSBG_P254_sin1_gap_n1_c1570_BW3T.txt', 'PSBG_P254_sin1_gap_n2_c1570_BW3T.txt'],
                     ['PSBG_P254_sin1_gap_n1_c1580_BW3T.txt', 'PSBG_P254_sin1_gap_n2_c1580_BW3T.txt'],
                     ['PSBG_P254_sin1_gap_n1_c1590_BW3T.txt', 'PSBG_P254_sin1_gap_n2_c1590_BW3T.txt'],
                     ['PSBG_P254_sin1_gap_n1_c1600_BW3T.txt', 'PSBG_P254_sin1_gap_n2_c1600_BW3T.txt'],
                     ['PSBG_P254_sin1_gap_n1_c1610_BW3T.txt', 'PSBG_P254_sin1_gap_n2_c1610_BW3T.txt'],
                     ['PSBG_P254_sin1_gap_n1_c1580_BW6T.txt', 'PSBG_P254_sin1_gap_n2_c1580_BW6T.txt'],
                     ['PSBG_P254_sin1_gap_n1_c1590_BW6T.txt', 'PSBG_P254_sin1_gap_n2_c1590_BW6T.txt'],
                     ['PSBG_P254_sin1_gap_n1_c1600_BW6T.txt', 'PSBG_P254_sin1_gap_n2_c1600_BW6T.txt'],
                     ['PSBG_P254_sin1_gap_n1_c1610_BW6T.txt', 'PSBG_P254_sin1_gap_n2_c1610_BW6T.txt']]

psbg_name_list_P254 = ['PSBG_c1570_BW3T_P254',
             'PSBG_c1580_BW3T_P254',
             'PSBG_c1590_BW3T_P254',
             'PSBG_c1600_BW3T_P254',
             'PSBG_c1610_BW3T_P254',
             'PSBG_c1580_BW6T_P254',
             'PSBG_c1590_BW6T_P254',
             'PSBG_c1600_BW6T_P254',
             'PSBG_c1610_BW6T_P254']

# psbg parameters
# first psbg cell
cell_x_P260 = cell_x_P254#int
cell_y_P260 = cell_y - group_space_y

param_psbg_P260 =  {"N_period": 1000,
                "SWG_period": 0.260, # micron
                "wg_width": 0.4,
                "Wing1_width": 0.130,
                "Wing2_width": 0.130,
                "leng_swg": 0.130}

filename_list_P260 = [['PSBG_P260_sin1_gap_n1_c1635_BW3T.txt', 'PSBG_P260_sin1_gap_n2_c1635_BW3T.txt'],
                     ['PSBG_P260_sin1_gap_n1_c1630_BW3T.txt', 'PSBG_P260_sin1_gap_n2_c1630_BW3T.txt'],
                     ['PSBG_P260_sin1_gap_n1_c1620_BW3T.txt', 'PSBG_P260_sin1_gap_n2_c1620_BW3T.txt'],
                     ['PSBG_P260_sin1_gap_n1_c1615_BW3T.txt', 'PSBG_P260_sin1_gap_n2_c1615_BW3T.txt'],
                     ['PSBG_P260_sin1_gap_n1_c1610_BW3T.txt', 'PSBG_P260_sin1_gap_n2_c1610_BW3T.txt'],
                     ['PSBG_P260_sin1_gap_n1_c1600_BW3T.txt', 'PSBG_P260_sin1_gap_n2_c1600_BW3T.txt'],
                     ['PSBG_P260_sin1_gap_n1_c1610_BW6T.txt', 'PSBG_P260_sin1_gap_n2_c1610_BW6T.txt'],
                     ['PSBG_P260_sin1_gap_n1_c1620_BW6T.txt', 'PSBG_P260_sin1_gap_n2_c1620_BW6T.txt'],
                     ['PSBG_P260_sin1_gap_n1_c1630_BW6T.txt', 'PSBG_P260_sin1_gap_n2_c1630_BW6T.txt']]

psbg_name_list_P260 = ['PSBG_c1635_BW3T_P260',
             'PSBG_c1630_BW3T_P260',
             'PSBG_c1620_BW3T_P260',
             'PSBG_c1615_BW3T_P260',
             'PSBG_c1610_BW3T_P260',
             'PSBG_c1600_BW3T_P260',
             'PSBG_c1610_BW6T_P260',
             'PSBG_c1620_BW6T_P260',
             'PSBG_c1630_BW6T_P260']

################# CBG #########################
dir_path_cbg = 'F:\\SWG\\Code\\CBG_SWG\\layout_list\\gap_list_Apr6_2021_refab\\'

n_x_cbg = 1
n_y_cbg = 4
# cbg p240
cell_x_P240_CBG = cell_x #int
cell_y_P240_CBG = cell_y - 3*group_space_y

param_cbg_P240 =  {"N_period": 20000,
                "SWG_period": 0.240, # micron
                "wg_width": 0.4,
                "Wing1_width": 0.120,
                "Wing2_width": 0.120,
                "leng_swg": 0.120}

filename_list_cbg_P240 = [['CBG_P240_gap_n1_v3gs2_c1520_BW_10000.txt', 'CBG_P240_gap_n2_v3gs2_c1520_BW_10000.txt'],
             ['CBG_P240_gap_n1_v4gs2_c1520_BW_10000.txt', 'CBG_P240_gap_n2_v4gs2_c1520_BW_10000.txt'],
             ['CBG_P240_gap_n1_v4gs3_c1520_BW_10000.txt', 'CBG_P240_gap_n2_v4gs3_c1520_BW_10000.txt'],
             ['CBG_P240_gap_n1_v4noapd_c1520_BW_10000.txt', 'CBG_P240_gap_n1_v4noapd_c1520_BW_10000.txt']]

cbg_name_list_P240 = ['CBG_c1520_v3gs2_P240_N10000',
                      'CBG_c1520_v4gs2_P240_N10000',
                      'CBG_c1520_v4gs3_P240_N10000',
                      'CBG_c1520_v4noapd_P240_N10000']

# cbg p245
cell_x_P245_CBG = cell_x #int
cell_y_P245_CBG = cell_y_P240_CBG - n_y_cbg*yspace

param_cbg_P245 =  {"N_period": 20000,
                "SWG_period": 0.245, # micron
                "wg_width": 0.4,
                "Wing1_width": 0.122,
                "Wing2_width": 0.122,
                "leng_swg": 0.122}

filename_list_cbg_P245 = [['CBG_P245_gap_n1_v3gs2_c1540_BW_10000.txt', 'CBG_P245_gap_n2_v3gs2_c1540_BW_10000.txt'],
             ['CBG_P245_gap_n1_v4gs2_c1540_BW_10000.txt', 'CBG_P245_gap_n2_v4gs2_c1540_BW_10000.txt'],
             ['CBG_P245_gap_n1_v4gs3_c1540_BW_10000.txt', 'CBG_P245_gap_n2_v4gs3_c1540_BW_10000.txt'],
             ['CBG_P245_gap_n1_v4noapd_c1540_BW_10000.txt', 'CBG_P245_gap_n1_v4noapd_c1540_BW_10000.txt']]

cbg_name_list_P245 = ['CBG_c1540_v3gs2_P245_N10000',
                      'CBG_c1540_v4gs2_P245_N10000',
                      'CBG_c1540_v4gs3_P245_N10000',
                      'CBG_c1540_v4noapd_P245_N10000']

# cbg p250
cell_x_P250_CBG = cell_x #int
cell_y_P250_CBG = cell_y_P245_CBG - n_y_cbg*yspace

param_cbg_P250 =  {"N_period": 20000,
                "SWG_period": 0.250, # micron
                "wg_width": 0.4,
                "Wing1_width": 0.125,
                "Wing2_width": 0.125,
                "leng_swg": 0.125}

filename_list_cbg_P250 = [['CBG_P250_gap_n1_v3gs2_c1565_BW_10000.txt', 'CBG_P250_gap_n2_v3gs2_c1565_BW_10000.txt'],
                        ['CBG_P250_gap_n1_v4gs2_c1565_BW_10000.txt', 'CBG_P250_gap_n2_v4gs2_c1565_BW_10000.txt'],
                        ['CBG_P250_gap_n1_v4gs3_c1565_BW_10000.txt', 'CBG_P250_gap_n2_v4gs3_c1565_BW_10000.txt'],
                        ['CBG_P250_gap_n1_v4noapd_c1565_BW_10000.txt', 'CBG_P250_gap_n1_v4noapd_c1565_BW_10000.txt']]

cbg_name_list_P250 = ['CBG_c1565_v3gs2_P250_N10000',
                      'CBG_c1565_v4gs2_P250_N10000',
                      'CBG_c1565_v4gs3_P250_N10000',
                      'CBG_c1565_v4noapd_P250_N10000']

# cbg p254
cell_x_P254_CBG = cell_x #int
cell_y_P254_CBG = cell_y_P250_CBG - n_y_cbg*yspace

param_cbg_P254 =  {"N_period": 20000,
                "SWG_period": 0.254, # micron
                "wg_width": 0.4,
                "Wing1_width": 0.127,
                "Wing2_width": 0.127,
                "leng_swg": 0.127}

filename_list_cbg_P254 = [['CBG_P254_gap_n1_v3gs2_c1585_BW_10000.txt', 'CBG_P254_gap_n2_v3gs2_c1585_BW_10000.txt'],
                        ['CBG_P254_gap_n1_v4gs2_c1585_BW_10000.txt', 'CBG_P254_gap_n2_v4gs2_c1585_BW_10000.txt'],
                        ['CBG_P254_gap_n1_v4gs3_c1585_BW_10000.txt', 'CBG_P254_gap_n2_v4gs3_c1585_BW_10000.txt'],
                        ['CBG_P254_gap_n1_v4noapd_c1585_BW_10000.txt', 'CBG_P254_gap_n1_v4noapd_c1585_BW_10000.txt']]

cbg_name_list_P254 = ['CBG_c1585_v3gs2_P254_N10000',
                      'CBG_c1585_v4gs2_P254_N10000',
                      'CBG_c1585_v4gs3_P254_N10000',
                      'CBG_c1585_v4noapd_P254_N10000']

# cbg p260
cell_x_P260_CBG = cell_x #int
cell_y_P260_CBG = cell_y_P254_CBG - n_y_cbg*yspace

param_cbg_P260 =  {"N_period": 20000,
                "SWG_period": 0.260, # micron
                "wg_width": 0.4,
                "Wing1_width": 0.130,
                "Wing2_width": 0.130,
                "leng_swg": 0.130}

filename_list_cbg_P260 = [['CBG_P260_gap_n1_v4noapd_c1615_BW_10000.txt', 'CBG_P260_gap_n2_v4noapd_c1615_BW_10000.txt'],
                        ['CBG_P260_gap_n1_v4gs3_c1615_BW_10000.txt', 'CBG_P260_gap_n2_v4gs3_c1615_BW_10000.txt'],
                        ['CBG_P260_gap_n1_v4gs2_c1615_BW_10000.txt', 'CBG_P260_gap_n2_v4gs2_c1615_BW_10000.txt'],
                        ['CBG_P260_gap_n1_v3gs2_c1615_BW_10000.txt', 'CBG_P260_gap_n2_v3gs2_c1615_BW_10000.txt']]

cbg_name_list_P260 = ['CBG_c1615_v4noapd_P260_N10000',
                      'CBG_c1615_v4gs3_P260_N10000',
                      'CBG_c1615_v4gs2_P260_N10000',
                      'CBG_c1615_v3gs2_P260_N10000'
                      ]

############## tilted SWG BG param #################
tswgbg_name = 'test'
n_x_tswgbg = 6
n_y_tswgbg = 11

cell_x_tswgbg = cell_x_P260 + 2*xspace #int
cell_y_tswgbg= cell_y_P260 - yspace

param_tswgbg = {"N_period": [1000],
                 "SWG_period": [0.258, 0.256],
                 "wg_width": [0.485],
                 "duty": [0.5],
                 "theta": [0,15,30,40,41,42,43,44,45,46,47,48,49,50,52],
                 "Wing1_width": [0.1],
                 "Wing2_width": [0.1],
                 "distance_2": [0.235,0.245]}

vgc_mode = 'TM'

# draw layout
draft = DrawLayout()
# draw FP box
tc = draft.FP(x_span, y_span, x_c, y_c)
# draw PSBG arrays
tc = draft.SWG_PSBG_array(psbg_name_list_P240, dir_path, filename_list_P240, param_psbg_P240, cell_x_P240, cell_y_P240, n_x, n_y, xspace, yspace)
tc = draft.SWG_PSBG_array(psbg_name_list, dir_path, filename_list, param_psbg1, cell_x_P245, cell_y_P245, n_x, n_y, xspace, yspace)
tc = draft.SWG_PSBG_array(psbg_name_list_P250, dir_path, filename_list_P250, param_psbg_P250, cell_x_P250, cell_y_P250, n_x, n_y, xspace, yspace)
tc = draft.SWG_PSBG_array(psbg_name_list_P254, dir_path, filename_list_P254, param_psbg_P254, cell_x_P254, cell_y_P254, n_x+1, n_y, xspace, yspace)
tc = draft.SWG_PSBG_array(psbg_name_list_P260, dir_path, filename_list_P260, param_psbg_P260, cell_x_P260, cell_y_P260, n_x+1, n_y, xspace, yspace)
# draw CBG arrays
tc = draft.SWG_PSBG_array(cbg_name_list_P240, dir_path_cbg, filename_list_cbg_P240, param_cbg_P240, cell_x_P240_CBG, cell_y_P240_CBG, n_x_cbg, n_y_cbg, xspace, yspace)
tc = draft.SWG_PSBG_array(cbg_name_list_P245, dir_path_cbg, filename_list_cbg_P245, param_cbg_P245, cell_x_P245_CBG, cell_y_P245_CBG, n_x_cbg, n_y_cbg, xspace, yspace)
tc = draft.SWG_PSBG_array(cbg_name_list_P250, dir_path_cbg, filename_list_cbg_P250, param_cbg_P250, cell_x_P250_CBG, cell_y_P250_CBG, n_x_cbg, n_y_cbg, xspace, yspace)
tc = draft.SWG_PSBG_array(cbg_name_list_P254, dir_path_cbg, filename_list_cbg_P254, param_cbg_P254, cell_x_P254_CBG, cell_y_P254_CBG, n_x_cbg, n_y_cbg, xspace, yspace)
tc = draft.SWG_PSBG_array(cbg_name_list_P260, dir_path_cbg, filename_list_cbg_P260, param_cbg_P260, cell_x_P260_CBG, cell_y_P260_CBG, n_x_cbg, n_y_cbg, xspace, yspace)
# draw tilted SWG BGs

tc = draft.Tilted_SWGBG_TETM_array(param_tswgbg, cell_x_tswgbg, cell_y_tswgbg, n_x_tswgbg, n_y_tswgbg, xspace, yspace)
tc.write("F:\Fab\ANT\Apr_6_2021_refab\\Apr_6_2021_refab_psbg_cbg_tiltbg.gds")
