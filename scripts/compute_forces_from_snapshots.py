import numpy as np
import readsnapHDF5 as rs

time_offset=0


if __name__ == '__main__':

    #binary and disk properties
    qb=float(sys.argv[1])
    eb=float(sys.argv[2])
    h=0.1
    alpha=0.1

    print eb

    if (len(sys.argv) < 4): init_snap = 0
    else:   init_snap = int(sys.argv[3]) 
    if (len(sys.argv) < 5): final_snap = 1001
    else:  final_snap = int(sys.argv[4]) 
    if (len(sys.argv) < 6): snap_step = 1
    else:  snap_step = int(sys.argv[5]) 
    
    snap_list = np.arange(init_snap,final_snap+1,snap_step)

    # paths to files
    run_path = "/data2/djmunoz/CIRCUMBINARY_DISKS_2D/RESTART_PLUTO_RUNS/"
    
    run_base ="restart-3000-pluto-woboundary-standardres-binary"

    run_name= run_base+("_q%.1f_e%.1f_h%.1f_alpha%.1f/" % (qb,eb,h,alpha))

    orbit_range = []
    for snap in [init_snap,final_snap]:
        filename=directory+base+snap_base+str(snap).zfill(3)
        header = rs.snapshot_header(filename)
        time = header.time + time_offset
        orbit_range.append(int(time/(2*np.pi)))
        acc = rs.read_block
