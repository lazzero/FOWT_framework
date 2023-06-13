import gmsh
import sys
import math
import numpy as np

def create_hydraspar_mesh_notclipped(name_file, draft, interface_TwSL,
                         spar_central1_radius, spar_central1_height,
                         spar_central2_radius, spar_central2_height,
                         heave_plate_radius, heave_plate_thickness,
                         spar_lateral_radius, spar_lateral_length,
                         spar_lateral_N, spar_lateral_angle,
                         alfa_init,
                         mesh_size_min, mesh_size_max,
                         *,show):

    gmsh.initialize()
    gmsh.option.setNumber("Mesh.MshFileVersion", 2)
    gmsh.model.add("temporary")

    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size_min)
    gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size_max)

    # define first cylinder
    XC1 = 0
    YC1 = 0
    ZC1 = interface_TwSL
    dxC1 = 0
    dyC1 = 0
    dzC1 = -spar_central1_height

    # define cone junction
    XCo1 = 0
    YCo1 = 0
    ZCo1 = ZC1 + dzC1
    dxCo1 = 0
    dyCo1 = 0
    dzCo1 = -((draft+interface_TwSL)-(spar_central1_height+spar_central2_height+heave_plate_thickness))

    # define second cylinder

    XC2 = 0
    YC2 = 0
    ZC2 = ZCo1 + dzCo1
    dxC2 = 0
    dyC2 = 0
    dzC2 = -spar_central2_height

    # define heave plate

    XHP1 = 0
    YHP1 = 0
    ZHP1 = ZC2+dzC2
    dxHP1 = 0
    dyHP1 = 0
    dzHP1 = -heave_plate_thickness

    RC1 = spar_central1_radius
    RC2 = spar_central2_radius
    RC3 = heave_plate_radius

    #define lateral cyilinder

    N = spar_lateral_N
    beta = spar_lateral_angle
    l = spar_lateral_length
    alfai = 360/N
    alfa_vec = alfa_init+np.linspace(alfai,360,N)

    #gmsh.model.occ.addBox(-R, -R, -R, 2 * R, 2 * R, 2 * R, 1)
    #gmsh.model.occ.addSphere(0, 0, 0, Rt, 2)
    C1 = gmsh.model.occ.addCylinder(XC1, YC1, ZC1, dxC1, dyC1, dzC1, RC1, tag = 1)
    Co1 = gmsh.model.occ.addCone(XCo1, YCo1, ZCo1, dxCo1,dyCo1,dzCo1,RC1,RC2, tag = 2)
    C2 = gmsh.model.occ.addCylinder(XC2, YC2, ZC2, dxC2, dyC2, dzC2, RC2, tag = 3)
    HP1 = gmsh.model.occ.addCylinder(XHP1, YHP1, ZHP1, dxHP1, dyHP1, dzHP1, RC3, tag = 4)

    Cl = np.empty(N)

    for i in range(0,N):

        alfa = alfa_vec[i]
        XCl = l*(0.07) *math.sin(math.radians(beta))*math.cos(math.radians(alfa))
        YCl = l*(0.07) *math.sin(math.radians(beta))*math.sin(math.radians(alfa))
        ZCl = -draft + l*(0.07)*math.cos(math.radians(beta))
        RCl = spar_lateral_radius
        dxCl = l*(0.93) *math.sin(math.radians(beta))*math.cos(math.radians(alfa))
        dyCl = l*(0.93) *math.sin(math.radians(beta))*math.sin(math.radians(alfa))
        dzCl = l*(0.93) *math.cos(math.radians(beta))

        Cl = np.append(Cl, [gmsh.model.occ.addCylinder(XCl,YCl,ZCl,dxCl,dyCl,dzCl,RCl, tag = 100+i)], axis=0)
        gmsh.model.occ.cut([(3,100+i)],[(3,3)],removeTool=False)
        gmsh.model.occ.remove([(3,100+i)], recursive=False)
        surfaces = gmsh.model.occ.getEntities(2)
        end = len(surfaces)
        print(surfaces)
        # if alfa == 180.0:
        #     gmsh.model.occ.remove([surfaces[(end-3)],surfaces[(end-2)],surfaces[(end-1)]], recursive=True)
        # else:
        #     # gmsh.model.occ.remove([surfaces[(end-2)],surfaces[(end-1)]], recursive=True)
        print('last is lateral spar lateral surface --------')
        print(surfaces)
        print('lines --------')
        lines = gmsh.model.occ.getEntities(1)
        print(lines)
        print('---------')

    # Create HP top surface
    gmsh.model.occ.cut([(2, 12)], [(2, 8)], 100)

    # Delete volumes and internal surfaces
    gmsh.model.occ.remove([(3,C1),(3,Co1),(3,C2),(3,HP1)], recursive=False)
    gmsh.model.occ.remove([(2,2),(2,3),(2,5),(2,6),(2,8),(2,9),(2,12)], recursive=True)
    print('lines @ this point --------')
    lines = gmsh.model.occ.getEntities(1)
    print(lines)
    print('---------')

    # create tower base
    gmsh.model.occ.addCurveLoop([3], 600)
    gmsh.model.occ.addPlaneSurface([600], 600)

    # # # Create sea level cut plane
    # # gmsh.model.occ.addPoint(20, 20, 0, tag = 100)
    # # gmsh.model.occ.addPoint(-20, 20, 0, tag = 200)
    # # gmsh.model.occ.addPoint(-20, -20, 0, tag = 300)
    # # gmsh.model.occ.addPoint(20, -20, 0, tag = 400)
    # # gmsh.model.occ.addLine(100, 200, tag = 100)
    # # gmsh.model.occ.addLine(300, 200, tag = 200)
    # # gmsh.model.occ.addLine(300, 400, tag = 300)
    # # gmsh.model.occ.addLine(400, 100, tag = 400)
    # # gmsh.model.occ.addCurveLoop([400, 100, -200, 300], 200)
    # # gmsh.model.occ.addPlaneSurface([200], 200)


    # # p = gmsh.model.occ.cut([(2, 1)], [(2, 200)])
    # # gmsh.model.occ.remove([p[0][2]],recursive=True)

    # # for i in range(0,N):
    # #     t = gmsh.model.occ.cut([(2, 13+i)], [(2, 200)])
    # #     gmsh.model.occ.remove([t[0][2]],recursive=True)

    # # # remove sea level
    # # gmsh.model.occ.remove([(2,200)],recursive=True)

    # print('lines @ end --------')
    # lines = gmsh.model.occ.getEntities(1)
    # print(lines)
    # print('---------')

    # # holes = np.empty(N)
    # # cut_surf = []
    # # k=0
    # # for i in range(0,N):
    # #     alfa = alfa_vec[i]
    # #     if alfa == 180.0:
    # #         l1 = gmsh.model.occ.cut([(1, 8)], [(1,404+k+3*i)],removeTool=False)
    # #         a = gmsh.model.occ.addCurveLoop([404+k+3*i,l1[0][1][1]])
    # #         b = gmsh.model.occ.addCurveLoop([405+k+3*i,l1[0][1][1]])
    # #         f1 = gmsh.model.occ.addPlaneSurface([a])
    # #         f2 = gmsh.model.occ.addPlaneSurface([b])
    # #         print('f1---------')
    # #         print(f1)
    # #         print('f2---------')
    # #         print(f2)
    # #         holes[i] = f1
    # #         c = gmsh.model.occ.cut([(2, 7)], [(2,f1),(2,f2)],removeObject=False)
    # #         print('cut_surf---------')
    # #         print(c)
    # #         print('---------')
    # #         cut_surf.append(c[0])
    # #         k=1
    # #     else:
    # #         a = gmsh.model.occ.addCurveLoop([404+k+3*i])
    # #         f = gmsh.model.occ.addPlaneSurface([a])
    # #         print('f---------')
    # #         print(f)
    # #         holes[i] = f
    # #         c = gmsh.model.occ.cut([(2, 7)], [(2,f)],removeObject=False)
    # #         print('cut_surf---------')
    # #         print(c)
    # #         print('---------')
    # #         cut_surf.append(c[0])



    # # print('holes---------')
    # # print(holes)
    # # print('---------')


    # # gmsh.model.occ.remove([(2,7)],recursive=True)

    # # intersect_surf = []
    # # intersect_surf.append(cut_surf[0])

    # # for i in range(0,N-1):
    # #     intersect_surf.append(gmsh.model.occ.intersect(intersect_surf[i],cut_surf[i+1])[0])

    # create mesh
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.recombine()
    gmsh.model.mesh.recombine()

    # Create name of file
    ext_file = ".msh"
    complete_name = "".join([name_file,ext_file])
    # ... and save it to disk
    gmsh.write(complete_name)

    # show mesh if needed
    if '-nopopup' not in sys.argv:
        if show:
            gmsh.fltk.run()

    gmsh.finalize()

    return True

# %% with braaces
def create_hydraspar_brace_mesh_notclipped(name_file, draft, interface_TwSL,
                         spar_central1_radius, spar_central1_height,
                         spar_central2_radius, spar_central2_height,
                         heave_plate_radius, heave_plate_thickness,
                         spar_lateral_radius, spar_lateral_length,
                         spar_lateral_N, spar_lateral_angle,
                         brace_height, brace_radius, brace_angle,
                         alfa_init,
                         mesh_size_min, mesh_size_max,
                         *,show):

    gmsh.initialize()
    gmsh.option.setNumber("Mesh.MshFileVersion", 2)
    gmsh.model.add("temporary")

    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size_min)
    gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size_max)

    # define first cylinder
    XC1 = 0
    YC1 = 0
    ZC1 = interface_TwSL
    dxC1 = 0
    dyC1 = 0
    dzC1 = -spar_central1_height

    # define cone junction
    XCo1 = 0
    YCo1 = 0
    ZCo1 = ZC1 + dzC1
    dxCo1 = 0
    dyCo1 = 0
    dzCo1 = -((draft+interface_TwSL)-(spar_central1_height+spar_central2_height+heave_plate_thickness))

    # define second cylinder

    XC2 = 0
    YC2 = 0
    ZC2 = ZCo1 + dzCo1
    dxC2 = 0
    dyC2 = 0
    dzC2 = -spar_central2_height

    # define heave plate

    XHP1 = 0
    YHP1 = 0
    ZHP1 = ZC2+dzC2
    dxHP1 = 0
    dyHP1 = 0
    dzHP1 = -heave_plate_thickness

    RC1 = spar_central1_radius
    RC2 = spar_central2_radius
    RC3 = heave_plate_radius

    #define lateral cyilinder

    N = spar_lateral_N
    beta = spar_lateral_angle
    l = spar_lateral_length
    alfai = 360/N
    alfa_vec = alfa_init+np.linspace(alfai,360,N)

    #gmsh.model.occ.addBox(-R, -R, -R, 2 * R, 2 * R, 2 * R, 1)
    #gmsh.model.occ.addSphere(0, 0, 0, Rt, 2)
    C1 = gmsh.model.occ.addCylinder(XC1, YC1, ZC1, dxC1, dyC1, dzC1, RC1, tag = 1)
    Co1 = gmsh.model.occ.addCone(XCo1, YCo1, ZCo1, dxCo1,dyCo1,dzCo1,RC1,RC2, tag = 2)
    C2 = gmsh.model.occ.addCylinder(XC2, YC2, ZC2, dxC2, dyC2, dzC2, RC2, tag = 3)
    HP1 = gmsh.model.occ.addCylinder(XHP1, YHP1, ZHP1, dxHP1, dyHP1, dzHP1, RC3, tag = 4)


    Cl = np.empty(N)
    Cb = np.empty(N)

    for i in range(0,N):

        alfa = alfa_vec[i]
        XCl = l*(0.07) *math.sin(math.radians(beta))*math.cos(math.radians(alfa))
        YCl = l*(0.07) *math.sin(math.radians(beta))*math.sin(math.radians(alfa))
        ZCl = -draft + l*(0.07)*math.cos(math.radians(beta))
        RCl = spar_lateral_radius
        dxCl = l*(0.93) *math.sin(math.radians(beta))*math.cos(math.radians(alfa))
        dyCl = l*(0.93) *math.sin(math.radians(beta))*math.sin(math.radians(alfa))
        dzCl = l*(0.93) *math.cos(math.radians(beta))
        Cl = np.append(Cl, [gmsh.model.occ.addCylinder(XCl,YCl,ZCl,dxCl,dyCl,dzCl,RCl, tag = 100+i)], axis=0)

        # Brace
        lb=brace_height*(math.tan(math.radians(brace_angle))*math.tan(math.radians(beta))/(math.tan(math.radians(brace_angle))+math.tan(math.radians(beta))))/math.sin(math.radians(brace_angle))
        XCb=0
        YCb=0
        ZCb=-draft+brace_height
        dxCb=lb*math.sin(math.radians(brace_angle))*math.cos(math.radians(alfa))
        dyCb=lb*math.sin(math.radians(brace_angle))*math.sin(math.radians(alfa))
        dzCb=-lb*math.cos(math.radians(brace_angle))
        Cb = np.append(Cb, [gmsh.model.occ.addCylinder(XCb, YCb, ZCb, dxCb, dyCb, dzCb, brace_radius, tag = 200+i)], axis=0)

        surfaces = gmsh.model.occ.getEntities(2)
        end = len(surfaces)
        gmsh.model.occ.cut([(3,100+i)],[(3,3)],removeTool=False)
        gmsh.model.occ.cut([(3,200+i)],[(3,1)],removeTool=False)
        ppp=gmsh.model.occ.cut([(3,200+i)],[(3,100+i)],removeTool=False)
        print('...... ppp =',ppp)
        gmsh.model.occ.remove([(3,100+i)], recursive=False)
        gmsh.model.occ.remove([(3,200+i)], recursive=False)

    # Create HP top surface
    gmsh.model.occ.cut([(2, 12)], [(2, 8)], 100)

    # Delete volumes and internal surfaces
    gmsh.model.occ.remove([(3,C1),(3,Co1),(3,C2),(3,HP1)], recursive=False)
    gmsh.model.occ.remove([(2,2),(2,5),(2,6),(2,8),(2,9),(2,12)], recursive=True)
    print('lines @ this point --------')
    lines = gmsh.model.occ.getEntities(1)
    print(lines)
    print('---------')

    # # Create sea level cut plane
    # gmsh.model.occ.addPoint(20, 20, 0, tag = 100)
    # gmsh.model.occ.addPoint(-20, 20, 0, tag = 200)
    # gmsh.model.occ.addPoint(-20, -20, 0, tag = 300)
    # gmsh.model.occ.addPoint(20, -20, 0, tag = 400)
    # gmsh.model.occ.addLine(100, 200, tag = 100)
    # gmsh.model.occ.addLine(300, 200, tag = 200)
    # gmsh.model.occ.addLine(300, 400, tag = 300)
    # gmsh.model.occ.addLine(400, 100, tag = 400)
    # gmsh.model.occ.addCurveLoop([400, 100, -200, 300], 200)
    # gmsh.model.occ.addPlaneSurface([200], 200)


    # p = gmsh.model.occ.cut([(2, 1)], [(2, 200)])
    # gmsh.model.occ.remove([p[0][2]],recursive=True)

    # for i in range(0,N):
    #     t = gmsh.model.occ.cut([(2, 13+i)], [(2, 200)])
    #     gmsh.model.occ.remove([t[0][2]],recursive=True)

    # # remove sea level
    # gmsh.model.occ.remove([(2,200)],recursive=True)

    # print('lines @ end --------')
    # lines = gmsh.model.occ.getEntities(1)
    # print(lines)
    # print('---------')

    # holes = np.empty(N)
    # cut_surf = []
    # k=0
    # for i in range(0,N):
    #     alfa = alfa_vec[i]
    #     if alfa == 180.0:
    #         l1 = gmsh.model.occ.cut([(1, 8)], [(1,404+k+3*i)],removeTool=False)
    #         a = gmsh.model.occ.addCurveLoop([404+k+3*i,l1[0][1][1]])
    #         b = gmsh.model.occ.addCurveLoop([405+k+3*i,l1[0][1][1]])
    #         f1 = gmsh.model.occ.addPlaneSurface([a])
    #         f2 = gmsh.model.occ.addPlaneSurface([b])
    #         print('f1---------')
    #         print(f1)
    #         print('f2---------')
    #         print(f2)
    #         holes[i] = f1
    #         c = gmsh.model.occ.cut([(2, 7)], [(2,f1),(2,f2)],removeObject=False)
    #         print('cut_surf---------')
    #         print(c)
    #         print('---------')
    #         cut_surf.append(c[0])
    #         k=1
    #     else:
    #         a = gmsh.model.occ.addCurveLoop([404+k+3*i])
    #         f = gmsh.model.occ.addPlaneSurface([a])
    #         print('f---------')
    #         print(f)
    #         holes[i] = f
    #         c = gmsh.model.occ.cut([(2, 7)], [(2,f)],removeObject=False)
    #         print('cut_surf---------')
    #         print(c)
    #         print('---------')
    #         cut_surf.append(c[0])



    # print('holes---------')
    # print(holes)
    # print('---------')


    # gmsh.model.occ.remove([(2,7)],recursive=True)

    # intersect_surf = []
    # intersect_surf.append(cut_surf[0])

    # for i in range(0,N-1):
    #     intersect_surf.append(gmsh.model.occ.intersect(intersect_surf[i],cut_surf[i+1])[0])

    # create mesh
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.recombine()
    # gmsh.model.mesh.recombine()

    # Create name of file
    ext_file = ".msh"
    complete_name = "".join([name_file,ext_file])
    # ... and save it to disk
    gmsh.write(complete_name)

    # show mesh if needed
    if '-nopopup' not in sys.argv:
        if show:
            gmsh.fltk.run()

    gmsh.finalize()

    return True

# %%
if __name__=="__main__":

    # create_hydraspar_mesh_notclipped(name_file='test', draft=9, interface_TwSL=6,
    #                      spar_central1_radius=1, spar_central1_height=7,
    #                      spar_central2_radius=2.5, spar_central2_height=5,
    #                      heave_plate_radius=5, heave_plate_thickness=0.06,
    #                      spar_lateral_radius=1, spar_lateral_length=23.33585740291,
    #                      spar_lateral_N=4, spar_lateral_angle=50,
    #                      alfa_init=0,
    #                      mesh_size_min=0.3, mesh_size_max=0.3,
    #                      show=True)

    create_hydraspar_brace_mesh_notclipped(name_file='test', draft=9, interface_TwSL=6,
                         spar_central1_radius=1, spar_central1_height=7,
                         spar_central2_radius=2.5, spar_central2_height=5,
                         heave_plate_radius=5, heave_plate_thickness=0.06,
                         spar_lateral_radius=1, spar_lateral_length=23.33585740291,
                         spar_lateral_N=5, spar_lateral_angle=50,
                         brace_height=12, brace_radius=0.5, brace_angle=80,
                         alfa_init=0,
                         mesh_size_min=0.3, mesh_size_max=0.3,
                         show=True)
