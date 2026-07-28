import sys 
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import ConnectionPatch
from matplotlib.widgets import RectangleSelector
from sympy import *
from tkinter import Tk, Frame
from tkinter import Button, Label

class GUI:
    def __init__(self):
        self.create_GUI()

    def animation(self):
        with open("C:\\Software Development\\6 DOF Sim\\Data\\GNC_Data2.csv") as csvFile:
            dictReader = csv.DictReader(csvFile)
            listofDicts = list(dictReader)
        m_x_lst,m_y_lst,m_z_lst = [],[],[]
        t_x_lst,t_y_lst,t_z_lst = [],[],[]
        for item in listofDicts[:-100]:
            m_x_lst.append(float(item['Mx']))
            m_y_lst.append(float(item['My']))
            m_z_lst.append(float(item['Mz']))
            t_x_lst.append(float(item['Tx']))
            t_y_lst.append(float(item['Ty']))
            t_z_lst.append(float(item['Tz']))
        tx_o,ty_o,tz_o = t_x_lst[0],t_y_lst[0],t_z_lst[0]
        intercept_animation(tx_o,ty_o,tz_o,m_x_lst,m_y_lst,m_z_lst,t_x_lst,t_y_lst,t_z_lst)

    def target_range_graph(self):
        with open("C:\\Software Development\\6 DOF Sim\\Data\\GNC_Data2.csv") as csvFile:
            dictReader = csv.DictReader(csvFile)
            listofDicts = list(dictReader)
        time_lst,range_lst = [],[]
        for item in listofDicts[:-1000]:
            time_lst.append(float(item['Time']))
            range_lst.append(float(item['Range']))
        missile_to_target_range_graph(time_lst,range_lst)

    def quaternion_graphs(self):
        with open("C:\\Software Development\\6 DOF Sim\\Data\\GNC_Data2.csv") as csvFile:
            dictReader = csv.DictReader(csvFile)
            listofDicts = list(dictReader)
        time_lst,q1_lst,q2_lst,q3_lst,q4_lst,q_lst = [],[],[],[],[],[]
        for item in listofDicts[:-1000]:
            time_lst.append(float(item['Time']))
            q1_lst.append(float(item['q1']))
            q2_lst.append(float(item['q2']))
            q3_lst.append(float(item['q3']))
            q4_lst.append(float(item['q4']))
            q_lst.append(float(item['q']))
        Quaternion_graphs(time_lst,q1_lst,q2_lst,q3_lst,q4_lst,q_lst)

    def propulsion_graphs(self):
        with open("C:\\Software Development\\6 DOF Sim\\Data\\GNC_Data2.csv") as csvFile:
            dictReader = csv.DictReader(csvFile)
            listofDicts = list(dictReader)
        time_lst,thrust_lst,fuel_mass_lst = [],[],[]
        for item in listofDicts[:-200]:
            time_lst.append(float(item['Time']))
            thrust_lst.append(float(item['Thrust']))
            fuel_mass_lst.append(float(item['fuel mass']))
        Propulsion_graphs(time_lst,thrust_lst,fuel_mass_lst)  

    def Monte_Carlo_graphs(self):
        with open("C:\\Software Development\\6 DOF SIM\\Data\\Monte_Carlo_Data.csv") as csvFile:
            dictReader = csv.DictReader(csvFile)
            listofDicts = list(dictReader)
        miss_d = []
        for item in listofDicts:
            miss_d.append(float(item['miss'])*3.821)
        CP_vs_MD(miss_d)

    def inertial_graphs(self):
        with open("C:\\Software Development\\6 DOF SIM\\Data\\Kalman_data.csv") as csvFile:
            dictReader = csv.DictReader(csvFile)
            listofDicts = list(dictReader)
        time_lst,tpx_lst,px_lst,meas_px_lst = [],[],[],[]
        tvx_lst,vx_lst = [],[]
        ax_lst = []
        for item in listofDicts[:-1000]:
            time_lst.append(float(item['time']))
            tpx_lst.append(float(item['tpx']))
            px_lst.append(float(item['px']))
            meas_px_lst.append(float(item['mpx']))
            tvx_lst.append(float(item['tvx']))
            vx_lst.append(float(item['vx']))
            ax_lst.append(float(item['mbax']))
        Inertial_Kalman(time_lst,tpx_lst,px_lst,meas_px_lst,tvx_lst,vx_lst,ax_lst)

    def terminal_graphs(self):
        with open("C:\\Software Development\\6 DOF SIM\\Data\\Target_Kalman_data.csv") as csvFile:
            dictReader = csv.DictReader(csvFile)
            listofDicts = list(dictReader)
        time_lst,tpx_lst,px_lst = [],[],[]
        tvx_lst,vx_lst = [],[]
        tax_lst,ax_lst = [],[]
        for item in listofDicts[:-1000]:
            time_lst.append(float(item['time']))
            tpx_lst.append(float(item['tpx']))
            px_lst.append(float(item['px']))
            tvx_lst.append(float(item['tvx']))
            vx_lst.append(float(item['vx']))
            tax_lst.append(float(item['tax']))
            ax_lst.append(float(item['ax']))
        Terminal_Kalman(time_lst,tpx_lst,px_lst,tvx_lst,vx_lst,tax_lst,ax_lst)

    def spsa(self):
        with open("C:\\Software Development\\6 DOF SIM\\Data\\SPSA_data.csv") as csvFile:
            dictReader = csv.DictReader(csvFile)
            listofDicts = list(dictReader)
        iter_lst,cost_lst = [],[]
        t0,t1,t2 = [],[],[]
        for item in listofDicts:
            iter_lst.append(float(item['iter']))
            cost_lst.append(float(item['cost']))
            t0.append(float(item['theta0']))
            t1.append(float(item['theta1']))
            t2.append(float(item['theta2']))
        SPSA_graph(iter_lst,cost_lst,t0,t1,t2)

    def aero_force_graphs(self):
        with open("C:\\Software Development\\6 DOF SIM\\Data\\Aero_data.csv") as csvFile:
            dictReader = csv.DictReader(csvFile)
            listofDicts = list(dictReader)
        time_lst,X_lst,Y_lst,Z_lst = [],[],[],[]
        for item in listofDicts[:-1000]:
            time_lst.append(float(item['time']))
            X_lst.append(float(item['X']))
            Y_lst.append(float(item['Y']))
            Z_lst.append(float(item['Z']))
        Aero_forces(time_lst,X_lst,Y_lst,Z_lst)

    def aero_moment_graphs(self):
        with open("C:\\Software Development\\6 DOF SIM\\Data\\Aero_data.csv") as csvFile:
            dictReader = csv.DictReader(csvFile)
            listofDicts = list(dictReader)
        time_lst,L_lst,M_lst,N_lst = [],[],[],[]
        for item in listofDicts[:-1000]:
            time_lst.append(float(item['time']))
            L_lst.append(float(item['L']))
            M_lst.append(float(item['M']))
            N_lst.append(float(item['N']))
        Aero_moments(time_lst,L_lst,M_lst,N_lst)

    def aero_force_coeff_graphs(self):
        with open("C:\\Software Development\\6 DOF SIM\\Data\\Aero_data.csv") as csvFile:
            dictReader = csv.DictReader(csvFile)
            listofDicts = list(dictReader)
        time_lst,CD_lst,CC_lst,CL_lst = [],[],[],[]
        for item in listofDicts[:-1000]:
            time_lst.append(float(item['time']))
            CD_lst.append(float(item['CD']))
            CC_lst.append(float(item['CC']))
            CL_lst.append(float(item['CL']))
        Aero_force_coefficiets(time_lst,CD_lst,CC_lst,CL_lst)

    def aero_moment_coeff_graphs(self):
        with open("C:\\Software Development\\6 DOF SIM\\Data\\Aero_data.csv") as csvFile:
            dictReader = csv.DictReader(csvFile)
            listofDicts = list(dictReader)
        time_lst,Cl_lst,Cm_lst,Cn_lst = [],[],[],[]
        for item in listofDicts[:-1000]:
            time_lst.append(float(item['time']))
            Cl_lst.append(float(item['Cl']))
            Cm_lst.append(float(item['Cm']))
            Cn_lst.append(float(item['Cn']))
        Aero_moment_coefficiets(time_lst,Cl_lst,Cm_lst,Cn_lst)

    def create_GUI(self):
        root = Tk()
        frame = Frame(root)
        root.title("APN Missile GNC System Simulation")
        frame.pack()

        title_label = Label(frame, text='APN Missile GNC System Simulation')
        title_label.grid(row=0,column=0)   

        space_label = Label(frame,text=" ")
        for i in range(5):
            space_label.grid(row=1,column=i) 

        ani_label = Label(frame,text="Engagement Visualization")
        ani_label.grid(row=2,column=0)
        run_ani_button = Button(frame,text="Run Animation",width=14,command=self.animation)
        run_ani_button.grid(row=2,column=1)

        range_graph_label = Label(frame,text="Target Range Graph")
        range_graph_label.grid(row=3,column=0)
        range_graph_button = Button(frame,text="Show Graph",width=14,command=self.target_range_graph)
        range_graph_button.grid(row=3,column=1)

        quat_label = Label(frame,text="Quaternion Vector")
        quat_label.grid(row=4,column=0)
        quat_button = Button(frame,text="Show Graphs",width=14,command=self.quaternion_graphs)
        quat_button.grid(row=4,column=1)

        prop_label = Label(frame,text="Propulsion Data")
        prop_label.grid(row=5,column=0)
        prop_button = Button(frame,text="Show Graphs",width=14,command=self.propulsion_graphs)
        prop_button.grid(row=5,column=1)

        monte_label = Label(frame,text="Monte Carlo Data")
        monte_label.grid(row=6,column=0)
        monte_button = Button(frame,text="Show Graphs",width=14,command=self.Monte_Carlo_graphs)
        monte_button.grid(row=6,column=1)

        init_label = Label(frame,text="Inertial Data")
        init_label.grid(row=7,column=0)
        init_button = Button(frame,text="Show Graphs",width=14,command=self.inertial_graphs)
        init_button.grid(row=7,column=1)

        term_label = Label(frame,text="Terminal Data")
        term_label.grid(row=8,column=0)
        term_button = Button(frame,text="Show Graphs",width=14,command=self.terminal_graphs)
        term_button.grid(row=8,column=1)

        opt_label = Label(frame,text="SPSA Data")
        opt_label.grid(row=9,column=0)
        opt_button = Button(frame,text="Show Graphs",width=14,command=self.spsa)
        opt_button.grid(row=9,column=1)

        force_label = Label(frame,text="Aero Forces")
        force_label.grid(row=10,column=0)
        force_button = Button(frame,text="Show Graphs",width=14,command=self.aero_force_graphs)
        force_button.grid(row=10,column=1)

        moment_label = Label(frame,text="Aero Moments")
        moment_label.grid(row=11,column=0)
        moment_button = Button(frame,text="Show Graphs",width=14,command=self.aero_moment_graphs)
        moment_button.grid(row=11,column=1)

        fcoe_label = Label(frame,text="Aero Force Coeff.")
        fcoe_label.grid(row=12,column=0)
        fcoe_button = Button(frame,text="Show Graphs",width=14,command=self.aero_force_coeff_graphs)
        fcoe_button.grid(row=12,column=1)

        momentcoe_label = Label(frame,text="Aero Moment Coeff.")
        momentcoe_label.grid(row=13,column=0)
        momentcoe_button = Button(frame,text="Show Graphs",width=14,command=self.aero_moment_coeff_graphs)
        momentcoe_button.grid(row=13,column=1)

        root.mainloop()

def missile_walk(m_x_lst,m_y_lst,m_z_lst):
    start_pos = np.array([m_x_lst[0],m_y_lst[0],m_z_lst[0]])
    m_walk = np.array([start_pos])
    for i in range(0,len(m_x_lst),300):
        walk = np.array([m_x_lst[i],m_y_lst[i],m_z_lst[i]])
        m_walk = np.append(m_walk,walk)
        m_walk = np.array_split(m_walk, len(m_walk)/3)
        m_walk = np.array(m_walk)
    return m_walk

def target_walk(tx_o,ty_o,tz_o,t_x_lst,t_y_lst,t_z_lst):
    start_pos = np.array([tx_o,ty_o,tz_o])
    t_walk = np.array([start_pos])
    for i in range(0,len(t_x_lst),300):
        walk = np.array([t_x_lst[i],t_y_lst[i],t_z_lst[i]])
        t_walk = np.append(t_walk,walk)
        t_walk = np.array_split(t_walk, len(t_walk)/3)
        t_walk = np.array(t_walk)
    return t_walk

def update_lines(num,walks,lines):
    for line, walk in zip(lines,walks):
        line.set_data(walk[:num, :2].T)
        line.set_3d_properties(walk[:num, 2])
    return lines

def intercept_animation(tx_o,ty_o,tz_o,m_x_lst,m_y_lst,m_z_lst,t_x_lst,t_y_lst,t_z_lst):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    m_walks = [missile_walk(m_x_lst,m_y_lst,m_z_lst)]
    m_lines = [ax.plot([],[],[])[0] for _ in m_walks]
    t_walks = [target_walk(tx_o,ty_o,tz_o,t_x_lst,t_y_lst,t_z_lst)]
    t_lines = [ax.plot([],[],[])[0] for _ in t_walks]
    ax.legend([m_lines[0],t_lines[0]],['APN missile flight path','target flight path'])
    ax.set(xlim3d=(-max(t_x_lst),max(t_x_lst)),xlabel='X')
    ax.set(ylim3d=(-max(t_y_lst),max(t_x_lst)),ylabel='Y')
    ax.set(zlim3d=(-max(t_z_lst),max(t_z_lst)),zlabel='Z')
    missile_ani = animation.FuncAnimation(fig,update_lines,len(m_x_lst),fargs=(m_walks,m_lines),interval=100)
    target_ani = animation.FuncAnimation(fig,update_lines,len(t_x_lst),fargs=(t_walks,t_lines),interval=100)
    plt.show()

# creates a crosshair that follows your mouse and displays the values on the graph
class Crosshair:
    # creates the characteristics of the cross hair following the points as well as which data it tracks
    def __init__(self,ax,line):
        self.ax = ax
        # creates horizontal and vertcial lines
        self.h_line = ax.axhline(color="black", lw=.08, linestyle="--")
        self.v_line = ax.axvline(color="black", lw=.08, linestyle="--")
        # gets x and y coordinates of the intersection
        self.x, self.y = line.get_data()
        # intitializes index of crosshair location
        self._last_index = None
        # sets characteristics and location of displayed text
        self.text = ax.text(.5, .8, '', transform = ax.transAxes)

    # determines if the crosshair needs to be drawn on the figure
    def cross_hair(self,visible):
        # sets need_to_draw as the opposite of visible (Boolean)
        need_to_draw = self.h_line.get_visible() != visible
        # sets the lines and text of the crosshair visible/not visible based on parameter visible
        self.h_line.set_visible(visible)
        self.v_line.set_visible(visible)
        self.text.set_visible(visible)
        return need_to_draw
    
    # determins if the crosshair needs to be updated
    def mouse_moved(self,event):
        # if mouse not in the range of the graph
        if not event.inaxes:
            # tells program to not draw the crosshair
            self._last_index = None
            need_to_draw = self.cross_hair(False)
            if need_to_draw:
                self.ax.figure.canvas.draw()
        # if mouse is in the graph
        else:
            # draws crosshair at the nearest x data point to the mouse
            self.cross_hair(True)
            x,y = event.xdata, event.ydata
            # finds the index in the x data array where the x value of the mouse would be put to maintain order then sets the index variable to 
            # either that index or the length of the x data array -1 based on which value is smaller
            index = min(np.searchsorted(self.x,x), len(self.x)-1)
            # doesn't update crosshair if x value hasn't changed
            if index == self._last_index:
                return
            # updates location of the crosshair 
            self._last_index = index
            x = self.x[index]
            y = self.y[index]
            # sets the x y coordiantes that will be displayed to the user
            self.h_line.set_ydata([y])
            self.v_line.set_xdata([x])
            # rounds displayed cooridnates to 2 decimal places
            self.text.set_text(f'time={x:1.2f}s, range={y:1.2f}m')
            self.ax.figure.canvas.draw()

# creates class to handle actions controlled by the keyboard/mouse
class Keyboard:
    # initializes which plot and line the key will effect
    def __init__(self,ax,line):
        self.ax = ax 
        self.line = line
        self.x, self.y = line.get_data()

    # hides/shows lines on graphs
    def toggle_lines(self, event):
        sys.stdout.flush()
        if event.key == 't':
            # toggles the line between visible and not visible
            visible = self.line.get_visible()
            self.line.set_visible(not visible)
            self.ax.figure.canvas.draw() 

    # grpahs the selected zoomed in region
    def zoom_in(self, event):
        # creates the figure for the new zoomed graph and plots the data
        fig_zoom, ax_zoom = plt.subplots()
        ax_zoom.set(title='Zoom window', xlabel='time (s)', ylabel='range')
        plt.plot(self.x,self.y)
        # does nothing if lmb isn't the buttom pressed
        if event.button != 1:
            return
        # sets x y coordinates to location of the mouse click
        x, y = event.xdata, event.ydata
        # sets the limits of x y axis to create the zoom
        ax_zoom.set_xlim(x-1, x+1)
        ax_zoom.set_ylim(y-30,y+30)
        # draws and shows figure to the user
        plt.axhline(.5, color = "black", linestyle = "--")
        plt.show()

# creates class for higlighting secitons of data
class DataSelector:
    # initializes attributes
    def __init__(self,fig,ax,z_ax,line):
        self.ax = ax
        self.fig = fig
        self.z_ax = z_ax
        self.x, self.y = line.get_data()

    # draws highlighted portion of main graph
    def live_zoom_window(self,c_event,r_event):
        # collects the x y cooridinates where the mouse was pressed and released
        x1, y1 = c_event.xdata, c_event.ydata
        x2, y2 = r_event.xdata, r_event.ydata
        # sets the new x y axis limits of the zoom graph to the coordinates of the highlighted box
        self.z_ax.set_xlim(x1,x2)
        self.z_ax.set_ylim(y1,y2)
        # draws the new graph
        self.fig.canvas.draw()

    # creates the highlighted rectangle
    def create_selector(self):
        selector = RectangleSelector(self.ax, self.live_zoom_window,
                                     button=1, minspanx=5, minspany=5,
                                     spancoords='pixels', interactive=True)
        return selector
    
# creates a class for making a connecting line between two graphs
class ConnectingLine:
    # initializes attributes
    def __init__(self,fig,ax_1,ax_2):
        self.fig = fig
        self.ax_1 = ax_1
        self.ax_2 = ax_2
    
    # creates a connecting line object and returns it
    def create_con_line(self):
        con_line = ConnectionPatch((0,0),(0,0),
                                   "data","data",
                                   axesA=self.ax_1,axesB=self.ax_2,
                                   color='black', linestyle='--')
        return con_line

def missile_to_target_range_graph(time,LOS_m_to_t_range_lst):
    # creates arrays to graph on x y axis
    LOS_m_to_t_range_lst = np.array(LOS_m_to_t_range_lst)
    time_axis = np.array(time)
    # generates and formats the graph
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.set_title("Missile to target range")
    # plots the data on the graph and draws a horizontal line at desginated value
    data_1,  = plt.plot(time_axis,LOS_m_to_t_range_lst)
    line_1 = plt.axhline(1, color = "black", linestyle = "--")
    # creates an object letting the data be affected by keyboard/mouse actions
    zoom_1 = Keyboard(ax1,data_1)
    plt.legend([line_1],['Intercept threshold'])
    plt.ylabel("range (m)")
    plt.tick_params('x')
    plt.xlabel("Time (s)")
    # creates a crosshair object that follows the specified data on the axis
    cursor_1 = Crosshair(ax1, data_1)
    fig.canvas.mpl_connect('motion_notify_event', cursor_1.mouse_moved)
    # displays graphs to user
    plt.show()

def Propulsion_graphs(time,thrust_lst,fuel_mass_lst):
    thrust_lst = np.array(thrust_lst)
    fuel_mass_lst = np.array(fuel_mass_lst)
    time_axis = np.array(time)
    fig = plt.figure()
    ax1 = plt.subplot(211)
    ax1.set_title("Thrust")
    data1, = plt.plot(time_axis,thrust_lst)
    plt.ylabel('Thrust')
    plt.tick_params('x',labelbottom=False)
    ax2 = plt.subplot(212, sharex = ax1)
    ax2.set_title('Fuel Mass')
    data2, = plt.plot(time_axis,fuel_mass_lst)
    plt.ylabel('Fuel Mass (kg)')
    plt.xlabel('Time (s)')
    plt.show()

def Quaternion_graphs(time,q1_lst,q2_lst,q3_lst,q4_lst,q_lst):
    time_axis = np.array(time)
    q1_lst = np.array(q1_lst)
    q2_lst = np.array(q2_lst)
    q3_lst = np.array(q3_lst)
    q4_lst = np.array(q4_lst)
    q_lst = np.array(q_lst)
    fig = plt.figure()
    ax1 = plt.subplot(511)
    ax1.set_title("Quaternion Vector During Flight")
    line_1, = plt.plot(time_axis,q1_lst)
    plt.ylabel("q1 (scalar)", fontsize=8)
    plt.tick_params('x',labelbottom=False)
    ax2 = plt.subplot(512, sharex= ax1)
    line_2, = plt.plot(time_axis,q2_lst)
    plt.ylabel("q2 (i)", fontsize=8)
    plt.tick_params('x',labelbottom=False)
    ax3 = plt.subplot(513, sharex = ax1)
    line_3, = plt.plot(time_axis,q3_lst)
    plt.ylabel("q3 (j)", fontsize=8)
    plt.tick_params('x',labelbottom=False)
    ax4 = plt.subplot(514, sharex = ax1)
    line_4, = plt.plot(time_axis,q4_lst)
    plt.ylabel("q4 (k)", fontsize=8)
    plt.tick_params('x',labelbottom=False)
    ax5 = plt.subplot(515,sharex = ax1)
    line_5, = plt.plot(time_axis,q_lst)
    plt.ylabel("q (magnitude)", fontsize=8)
    plt.xlabel("Time (s)")
    q1_zoom_fig,q1_zoom_ax = plt.subplots(figsize=(5,2))
    q1_zoom_ax.set(title="Highlighted q1",xlabel='Time (s)',ylabel='q1 (scalar)')
    plt.plot(time_axis,q1_lst)
    q1_selector = DataSelector(q1_zoom_fig,ax1,q1_zoom_ax,line_1)
    q1_highlighter = q1_selector.create_selector()
    q2_zoom_fig,q2_zoom_ax = plt.subplots(figsize=(5,2))
    q2_zoom_ax.set(title="Highlighted q2",xlabel='Time (s)',ylabel='q2 (i)')
    plt.plot(time_axis,q2_lst)
    q2_selector = DataSelector(q2_zoom_fig,ax2,q2_zoom_ax,line_2)
    q2_highlighter = q2_selector.create_selector()
    q3_zoom_fig,q3_zoom_ax = plt.subplots(figsize=(5,2))
    q3_zoom_ax.set(title="Highlighted q3",xlabel='Time (s)',ylabel='q3 (j)')
    plt.plot(time_axis,q3_lst)
    q3_selector = DataSelector(q3_zoom_fig,ax3,q3_zoom_ax,line_3)
    q3_highlighter = q3_selector.create_selector()
    q4_zoom_fig,q4_zoom_ax = plt.subplots(figsize=(5,2))
    q4_zoom_ax.set(title="Highlighted q4",xlabel="Time (s)",ylabel='q1 (scalar)')
    plt.plot(time_axis,q4_lst)
    q4_selector = DataSelector(q4_zoom_fig,ax4,q4_zoom_ax,line_4)
    q4_highlighter = q4_selector.create_selector()

    plt.show()

def CP_vs_MD(miss_d):
    miss_arr = np.array(miss_d)
    miss_arr = np.sort(miss_arr)

    N = len(miss_arr)
    cum_prob = np.arange(1, N+1) / N

    miss_arr = np.insert(miss_arr, 0, 0.0)
    cum_prob = np.insert(cum_prob, 0, 0.0)

    plt.plot(miss_arr, cum_prob)
    plt.xlabel('Miss Distance (ft)')
    plt.ylabel('Cumulative Probability')
    plt.title('Cumulative Probability vs Miss Distance')
    plt.grid(True)
    plt.show()

def Inertial_Kalman(time_lst,tpx_lst,px_lst,meas_px_lst,tvx_lst,vx_lst,ax_lst):
    time_arr = np.array(time_lst)
    tpx_arr = np.array(tpx_lst)
    px_arr = np.array(px_lst)
    mpx_arr = np.array(meas_px_lst)
    mask = ~np.isnan(mpx_arr)
    tvx_arr = np.array(tvx_lst)
    vx_arr = np.array(vx_lst)
    bax_arr = np.array(ax_lst)
    fig = plt.figure()
    ax1 = plt.subplot(311)
    ax1.set_title("Inertial Extended Kalman Filter")
    line1, = plt.plot(time_arr,tpx_arr)
    line2, = plt.plot(time_arr[mask],mpx_arr[mask],'r.-', alpha=0.7)
    line3, = plt.plot(time_arr,px_arr)
    #line2, = plt.scatter(time_arr,mpx_arr, s=10, color='red', alpha=0.6)
    plt.ylabel("x Postion")
    plt.legend([line1,line2,line3],["True Position","Measured position","Estimated Position"])
    plt.tick_params('x',labelbottom=False)
    ax2 = plt.subplot(312,sharex = ax1)
    vline1, = plt.plot(time_arr,tvx_arr)
    vline2, = plt.plot(time_arr,vx_arr)
    plt.ylabel("Velocity")
    plt.tick_params('x',labelbottom=False)
    plt.legend([vline1,vline2],["True Velocity","Estimated Velocity"])
    ax3 = plt.subplot(313,sharex = ax1)
    aline2, = plt.plot(time_arr,bax_arr,color = "red")
    plt.ylabel("Acceleration")
    plt.legend([aline2],["Measured Acceleration"])
    plt.xlabel("Time (s)")
    plt.show()

def Terminal_Kalman(time_lst,tpx_lst,px_lst,tvx_lst,vx_lst,tax_lst,ax_lst):
    time_arr = np.array(time_lst)
    tpx_arr = np.array(tpx_lst)
    px_arr = np.array(px_lst)
    tvx_arr = np.array(tvx_lst)
    vx_arr = np.array(vx_lst)
    tax_arr = np.array(tax_lst)
    ax_arr = np.array(ax_lst)
    fig = plt.figure()
    ax1 = plt.subplot(311)
    ax1.set_title("Terminal Extended Kalman Filter")
    line1, = plt.plot(time_arr,tpx_arr)
    line3, = plt.plot(time_arr,px_arr)
    #plt.scatter(time_arr,mpx_arr, s=10, color='red', alpha=0.6)
    plt.ylabel("x Postion")
    plt.legend([line1,line3],["True Position","Estimated Position"])
    plt.tick_params('x',labelbottom=False)
    ax2 = plt.subplot(312,sharex = ax1)
    vline1, = plt.plot(time_arr,tvx_arr)
    vline2, = plt.plot(time_arr,vx_arr)
    plt.ylabel("Velocity")
    plt.tick_params('x',labelbottom=False)
    plt.legend([vline1,vline2],["True Velocity","Estimated Velocity"])
    ax3 = plt.subplot(313,sharex = ax1)
    aline1, = plt.plot(time_arr,tax_arr)
    aline2, = plt.plot(time_arr,ax_arr)
    plt.ylabel("Acceleration")
    plt.legend([aline1,aline2],["True Acceleration","Estimated Acceleration"])
    plt.xlabel("Time (s)")
    plt.show()

def Aero_force_coefficiets(time_lst,CD_lst,CC_lst,CL_lst):
    time_arr = np.array(time_lst)
    CD_arr = np.array(CD_lst)
    CC_arr = np.array(CC_lst)
    CL_arr = np.array(CL_lst)
    fig = plt.figure()
    ax1 = plt.subplot(311)
    ax1.set_title("Aerodynamic Force Coefficients")
    line1, = plt.plot(time_arr,CD_arr)
    plt.ylabel("CD")
    plt.tick_params('x',labelbottom=False)
    ax2 = plt.subplot(312)
    line2, = plt.plot(time_arr,CC_arr)
    plt.ylabel("CC")
    plt.tick_params('x',labelbottom=False)
    ax3 = plt.subplot(313)
    line1, = plt.plot(time_arr,CL_arr)
    plt.ylabel("CL")
    plt.xlabel("Time (s)")
    plt.show()

def Aero_moment_coefficiets(time_lst,Cl_lst,Cm_lst,Cn_lst):
    time_arr = np.array(time_lst)
    Cl_arr = np.array(Cl_lst)
    Cm_arr = np.array(Cm_lst)
    Cn_arr = np.array(Cn_lst)
    fig = plt.figure()
    ax1 = plt.subplot(311)
    ax1.set_title("Aerodynamic Moment Coefficients")
    line1, = plt.plot(time_arr,Cl_arr)
    plt.ylabel("Cl")
    plt.tick_params('x',labelbottom=False)
    ax2 = plt.subplot(312)
    line2, = plt.plot(time_arr,Cm_arr)
    plt.ylabel("Cm")
    plt.tick_params('x',labelbottom=False)
    ax3 = plt.subplot(313)
    line1, = plt.plot(time_arr,Cn_arr)
    plt.ylabel("Cn")
    plt.xlabel("Time (s)")
    plt.show()

def Aero_forces(time_lst,X_lst,Y_lst,Z_lst):
    time_arr = np.array(time_lst)
    X_arr = np.array(X_lst)
    Y_arr = np.array(Y_lst)
    Z_arr = np.array(Z_lst)
    fig = plt.figure()
    ax1 = plt.subplot(311)
    ax1.set_title("Aerodynamic Forces")
    line1, = plt.plot(time_arr,X_arr)
    plt.ylabel("X")
    plt.tick_params('x',labelbottom=False)
    ax2 = plt.subplot(312)
    line2, = plt.plot(time_arr,Y_arr)
    plt.ylabel("Y")
    plt.tick_params('x',labelbottom=False)
    ax3 = plt.subplot(313)
    line1, = plt.plot(time_arr,Z_arr)
    plt.ylabel("Z")
    plt.xlabel("Time (s)")
    plt.show()

def Aero_moments(time_lst,L_lst,M_lst,N_lst):
    time_arr = np.array(time_lst)
    L_arr = np.array(L_lst)
    M_arr = np.array(M_lst)
    N_arr = np.array(N_lst)
    fig = plt.figure()
    ax1 = plt.subplot(311)
    ax1.set_title("Aerodynamic Moments")
    line1, = plt.plot(time_arr,L_arr)
    plt.ylabel("l")
    plt.tick_params('x',labelbottom=False)
    ax2 = plt.subplot(312)
    line2, = plt.plot(time_arr,M_arr)
    plt.ylabel("m")
    plt.tick_params('x',labelbottom=False)
    ax3 = plt.subplot(313)
    line1, = plt.plot(time_arr,N_arr)
    plt.ylabel("n")
    plt.xlabel("Time (s)")
    plt.show()

def SPSA_graph(iter_lst,cost_lst,t0,t1,t2):
    iter_arr = np.array(iter_lst)
    cost_arr = np.array(cost_lst)
    t0_arr = np.array(t0)
    t1_arr = np.array(t1)
    t2_arr = np.array(t2)
    fig = plt.figure()
    ax1 = plt.subplot(211)
    ax1.set_title("SPSA Parameter Optimization")
    plt.plot(iter_arr,cost_arr)
    plt.ylabel("Cost")
    plt.tick_params('x',labelbottom=False)
    ax2 = plt.subplot(212,sharex = ax1)
    line1, = plt.plot(iter_arr,t0_arr)
    line2, = plt.plot(iter_arr,t1_arr)
    line3, = plt.plot(iter_arr,t2_arr)
    plt.ylabel("Parameter Value")
    plt.xlabel("Iteration")
    plt.legend([line1,line2,line3],["Theta 0","Theta 1","Theta 2"])
    plt.show()


#engagement_simulation()

gui = GUI()