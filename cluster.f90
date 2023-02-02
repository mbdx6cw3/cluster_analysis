program cluster
! -----------------------------------------------------------------------------------
! This program performs cluster analysis on a molecular simulation trajectory.
! Clustering criterion is based purely on the distance between the centres of
! mass of defined molecules. It relies on a recursive subroutine (mol_cluster).
!
! Can be used to calculate:
!	- max number of clusters, mean and max cluster size and their time averages (data.dat)
!	- cluster size distributions (written to dist_out.dat)
!	- mean distance between clusters
!	- radial distribution functions (written to rdf_out.dat)
!
! Required input:
!	1) input.dat: contains input parameters (see example)
!	2) atoms.dat: a list of atom types and their masses 	
!	3) simulation trajector in .gro format. All molecules not involved
!	   in the definition of a cluster must be removed.
!
! To run the program, first compile global.f90, e.g.:
!	gfortran -c global.f90
! Then the main executable, e.g.:
!	gfortran -o cluster cluster.f90 global.f90
! 
! Enjoy!
!
! Christopher D. Williams (Dec. 2019).
! ----------------------------------------------------------------------------------
	use global
	implicit none

        ! --- read input file
        open(unit = 1, file = 'input.dat', status = 'old')
                read(unit = 1, nml = input)
        close(unit = 1)

	! --- allocate
	allocate(type(n_typ),mass(n_typ))

	! --- read atom info
        open(unit = 1, file = 'atoms.dat', status = 'old')
                read(unit = 1, nml = atoms)
        close(unit = 1)

        ! --- read in trajectory .gro file
        call read_trj

	! --- assign particle masses
        allocate(atm_mas(n_mol_atm))
	mol_mas = 0.0d0
	do i_atm = 1, n_mol_atm
		do i_typ = 1, n_typ
			if (atm_typ(i_atm).eq.type(i_typ)) then
				atm_mas(i_atm) = mass(i_typ)
                		mol_mas = mol_mas + atm_mas(i_atm)
			end if
		end do
	end do

        ! --- calculate number of frames to use cluster data
        n_frame = (end_frame - (start_frame - 1)) / space_frame
	allocate(x_mol(end_frame,n_mol),y_mol(end_frame,n_mol),z_mol(end_frame,n_mol),dist(end_frame,n_mol,n_mol))
	allocate(av_dist(end_frame))
        allocate(neighbour(end_frame,n_mol,max_neighbour),neigh_num(end_frame,n_mol),cluster_assign(end_frame,n_mol))
	allocate(cluster_size_dist(n_mol),cluster_time(end_frame,n_mol))
	cluster_assign (:,:) = .false.

        ! --- calculate distance between molecules centres-of-mass	
	call mol_com

	! --- calculate radial distribution function
	call calc_rdf

	! --- cluster size distribution analysis
	call cluster_dist

        ! --- cluster residence times
        call cluster_life

	! --- write output
	call write_out

end program

subroutine read_trj
        use global, only: n_mol_atm,n_mol,x,y,z,lx,ly,lz,end_frame,start_frame,space_frame,atm_typ
	implicit none

        ! --- variables
        character(len=13) :: A, nam
	character(len=2) :: t_atm
        integer :: atm_i, i_frame, i_mol, i_atm, n_atm

        ! --- read in initial coordinates for 3D periodic membrane model
        100 format (a13, a2, i5, f8.3, f8.3, f8.3)
        101 format (f10.5, f10.5, f10.5)

        ! --- allocate frame arrays
        allocate(lx(end_frame),ly(end_frame),lz(end_frame))

        open(unit=10,file='input.gro')
        do i_frame = 1, end_frame

                read(10,*) A
                read(10,*) n_atm

                ! --- allocate atom arrays in first frame
                if (i_frame.eq.1) then
		        n_mol_atm = n_atm / n_mol
                        allocate(x(end_frame,n_mol,n_mol_atm),y(end_frame,n_mol,n_mol_atm),z(end_frame,n_mol,n_mol_atm),&
			& atm_typ(n_mol_atm))
                end if

		! --- read in coordinates of first sheet
		do i_mol = 1, n_mol
			do i_atm = 1, n_mol_atm
				read(10,100) nam, t_atm, atm_i, x(i_frame, i_mol, i_atm), y(i_frame, i_mol, i_atm), z(i_frame, i_mol, i_atm)
				if ((i_frame.eq.1).and.(i_mol.eq.1)) then
					atm_typ(i_atm) = t_atm
                		end if
                	end do
		end do		
                read(10,*) lx(i_frame), ly(i_frame), lz(i_frame)

        end do
        close(unit=10)

end subroutine

subroutine mol_com
        use global, only: end_frame,start_frame,space_frame,n_mol,n_mol_atm,x,y,z,mol_mas,av_dist,&
	& atm_mas,x_mol,y_mol,z_mol,lx,ly,lz,dist,pi,dp,neigh_num,neighbour,cluster_cut
        implicit none

        ! --- variables
        integer :: i_frame, i_mol, j_mol, i_atm
        real(kind=dp) :: theta_x, theta_y, theta_z, theta_x_sum, theta_y_sum, theta_z_sum
        real(kind=dp) :: xi_x, xi_y, xi_z, xi_x_sum, xi_y_sum, xi_z_sum
        real(kind=dp) :: zeta_x, zeta_y, zeta_z, zeta_x_sum, zeta_y_sum, zeta_z_sum
	real(kind=dp) :: rx, ry, rz, sum_dist

        ! --- calculate com coordinates for each molecule
        do i_frame = start_frame, end_frame, space_frame

                ! --- loop over molecules
                do i_mol = 1, n_mol

                        xi_x_sum = 0.0d0
                        xi_y_sum = 0.0d0
                        xi_z_sum = 0.0d0

                        zeta_x_sum = 0.0d0
                        zeta_y_sum = 0.0d0
                        zeta_z_sum = 0.0d0

                        ! --- calculate sheet centre of mass position
                        ! --- coordinates must be treated as if they are on a sphere instead of a line
                        ! --- see wiki page "Center of mass" for details
                        do i_atm = 1, n_mol_atm

                                theta_x = 2.0d0 * pi * x(i_frame,i_mol,i_atm) / lx(i_frame)
                                xi_x = cos(theta_x)
                                zeta_x = sin(theta_x)
                                xi_x_sum = xi_x_sum + atm_mas(i_atm) * xi_x
                                zeta_x_sum = zeta_x_sum + atm_mas(i_atm) * zeta_x

                                theta_y = 2.0d0 * pi * y(i_frame,i_mol,i_atm) / ly(i_frame)
                                xi_y = cos(theta_y)
                                zeta_y = sin(theta_y)
                                xi_y_sum = xi_y_sum + atm_mas(i_atm) * xi_y
                                zeta_y_sum = zeta_y_sum + atm_mas(i_atm) * zeta_y

                                theta_z = 2.0d0 * pi * z(i_frame,i_mol,i_atm) / lz(i_frame)
                                xi_z = cos(theta_z)
                                zeta_z = sin(theta_z)
                                xi_z_sum = xi_z_sum + atm_mas(i_atm) * xi_z
                                zeta_z_sum = zeta_z_sum + atm_mas(i_atm) * zeta_z

                        end do

                        xi_x_sum = xi_x_sum / mol_mas
                        zeta_x_sum = zeta_x_sum / mol_mas
                        xi_y_sum = xi_y_sum / mol_mas
                        zeta_y_sum = zeta_y_sum / mol_mas
                        xi_z_sum = xi_z_sum / mol_mas
                        zeta_z_sum = zeta_z_sum / mol_mas

                        theta_x_sum = atan2(-1.0d0*zeta_x_sum,-1.0d0*xi_x_sum) + pi
                        x_mol(i_frame, i_mol) = (lx(i_frame) * theta_x_sum) / (2.0d0 * pi)

                        theta_y_sum = atan2(-1.0d0*zeta_y_sum,-1.0d0*xi_y_sum) + pi
                        y_mol(i_frame, i_mol) = (ly(i_frame) * theta_y_sum) / (2.0d0 * pi)

                        theta_z_sum = atan2(-1.0d0*zeta_z_sum,-1.0d0*xi_z_sum) + pi
                        z_mol(i_frame, i_mol) = (lz(i_frame) * theta_z_sum) / (2.0d0 * pi)

                end do

        end do

        ! --- calculate distance between all molecule coms
        do i_frame = start_frame, end_frame, space_frame

		sum_dist = 0.0d0
                ! --- loop over molecules
                do i_mol = 1, n_mol - 1

                        ! --- calculate distance to other molecules
                        do j_mol = i_mol + 1, n_mol

                                rx = x_mol(i_frame, i_mol) - x_mol(i_frame, j_mol)
                                ry = y_mol(i_frame, i_mol) - y_mol(i_frame, j_mol)
                                rz = z_mol(i_frame, i_mol) - z_mol(i_frame, j_mol)
                                rx = rx - lx(i_frame) * nint(rx/lx(i_frame))
                                ry = ry - ly(i_frame) * nint(ry/ly(i_frame))
                                rz = rz - lz(i_frame) * nint(rz/lz(i_frame))
                                dist(i_frame,i_mol,j_mol) = sqrt (rx * rx + ry * ry + rz * rz)
				sum_dist = sum_dist + dist(i_frame,i_mol,j_mol)
				if (dist(i_frame,i_mol,j_mol).lt.cluster_cut) then
                                        neigh_num(i_frame, i_mol) = neigh_num(i_frame, i_mol) + 1
                                        neigh_num(i_frame, j_mol) = neigh_num(i_frame, j_mol) + 1
                                        neighbour (i_frame, i_mol, neigh_num(i_frame, i_mol)) = j_mol
                                        neighbour (i_frame, j_mol, neigh_num(i_frame, j_mol)) = i_mol
				end if

                        end do

                end do
		av_dist(i_frame) = sum_dist / real((n_mol - 1) * (n_mol / 2))

        end do

end subroutine

subroutine calc_rdf
	use global, only: n_mol,x,y,z,lx,ly,lz,end_frame,space_frame,start_frame,n_frame,n_bin,dist,rdf_max,dp
        implicit none
	real(kind=dp), parameter :: pi = 3.14159265359d0
        real(kind=dp) :: bin_width, rx, ry, rz, bin_vol, max_bin, rcut
        integer :: i_frame, i_mol, j_mol, i_bin
	real(kind=dp), dimension (n_bin) :: gr
 
        200 format (f8.3, f11.3) 

	! --- RDF setup
        bin_width = rdf_max / real(n_bin)
	gr(:) = 0.0d0

        ! --- write histograms for Dirac delta function
        do i_frame = start_frame, end_frame, space_frame

		! --- calculate sheet COM - sheet COM RDFs
                do i_mol = 1, n_mol - 1

                        ! --- calculate distance to other molecules
                        do j_mol = i_mol + 1, n_mol

                                if (dist(i_frame,i_mol,j_mol).lt.rdf_max) then
                                	i_bin = ceiling(dist(i_frame,i_mol,j_mol)/bin_width)
					gr(i_bin) = gr(i_bin) + 1.0d0
				end if

			end do

		end do

        end do

	! --- calculate radial distribution function
	open(unit=20,file='rdf_out.dat',status='unknown')
        do i_bin = 1, n_bin
                bin_vol = (4.0d0/3.0d0) * pi * (i_bin**3 - (i_bin-1)**3) * bin_width**3
                gr(i_bin) = lx(1) * ly(1) * lz(1) * gr(i_bin) / bin_vol / real((n_mol-1)*(n_mol/2)) / real(n_frame)
                write(20,200) (i_bin - 0.5d0) * bin_width, gr(i_bin)
        end do
        close(unit=20)

end subroutine

subroutine cluster_dist
	use global, only: start_frame,end_frame,space_frame,n_frame,n_mol,dp,av_dist,&
			& cluster_count,cluster_size,cluster_assign,neigh_num,n_cluster,&
			& cluster_index
	implicit none
	integer :: size_count, i_mol, i_frame, n_mol_chk

        allocate(cluster_size(end_frame,n_mol),cluster_count(end_frame),cluster_index(end_frame,n_mol,n_mol))

        do i_frame = start_frame, end_frame, space_frame

		n_cluster = 0
        	n_mol_chk = 0
                do i_mol = 1, n_mol

                        ! --- if molecule has already been assigned to a cluster then go to next one
                        if (cluster_assign(i_frame, i_mol)) cycle

                        ! --- new size of this cluster
                        size_count = 1

                        ! --- update number of clusters for this frame
                        n_cluster = n_cluster + 1
                        
			! --- use recursive algorithm to determine size of this cluster from neighbour matrix
                        call mol_cluster(i_frame, i_mol, size_count)

                        ! --- size of this cluster
                        cluster_size(i_frame, n_cluster) = size_count
                        n_mol_chk = n_mol_chk + size_count

                end do

                ! --- check number of atoms
                if (n_mol_chk.ne.n_mol) then
                        print*, "ERROR - total number of atoms in cluster incorrect"
                        stop
                end if

                ! --- number of clusters in this frame
                cluster_count(i_frame) = n_cluster

	end do

end subroutine

subroutine cluster_life
	use global, only: start_frame,end_frame,space_frame,cluster_count,&
			& cluster_index,cluster_size,cluster_time,n_mol
	implicit none
	integer :: i_cluster, j_cluster, i_frame, i_mol, j_mol, i_time
	integer, allocatable, dimension(:,:) :: cluster_time_dist

	! --- allocate 
	allocate(cluster_time_dist(n_mol,end_frame))

	! --- initiate time counter
	cluster_time(start_frame,:) = 1

	! --- loop over all frames
	do i_frame = start_frame + 1, end_frame, space_frame

		! --- loop over all clusters in this frame
		do i_cluster = 1, cluster_count(i_frame)
	
                        ! --- reset timer
                        cluster_time(i_frame,i_cluster) = 1

			! --- loop over all clusters in previous frame
			do j_cluster = 1, cluster_count(i_frame-1)

				! --- if cluster is the same size of a cluster in previous frame it could be the same cluster
				if (cluster_size(i_frame,i_cluster).eq.cluster_size(i_frame-1,j_cluster)) then
					
					! --- if cluster has all the same molecules as cluster identified in previous frame then it is the same cluster
                                	do i_mol = 1, cluster_size(i_frame, i_cluster)

						do j_mol = 1, cluster_size(i_frame-1,j_cluster)

                                        		! --- if cluster is the same size check that molecules in it are the same
                                        		if (cluster_index(i_frame,i_cluster,i_mol).eq.&
							& cluster_index(i_frame-1,j_cluster,j_mol)) then

								! --- if final molecule in cluster matches then update cluster timer
								if (i_mol.eq.cluster_size(i_frame, i_cluster)) then

									! --- update current length of time this cluster has existed
                                                			cluster_time(i_frame,i_cluster) = cluster_time(i_frame-1,j_cluster) + 1

								end if
							
							end if
                                        	end do
					end do

				end if

			end do


			cluster_time_dist(cluster_size(i_frame,i_cluster),cluster_time(i_frame,i_cluster)) = &
                                        & cluster_time_dist(cluster_size(i_frame,i_cluster),cluster_time(i_frame,i_cluster)) + 1

		end do
		
	end do

	! --- write cluster time distributions
	do i_time = 1, 8000
		print*, i_time,cluster_time_dist(5,i_time) 
	end do

end subroutine

subroutine write_out
	use global, only: cluster_size,start_frame,end_frame,space_frame,cluster_count,&
			& av_dist,n_frame,n_mol,dp,cluster_index,out_style,cluster_size_dist,&
			& cluster_time,dt 
	implicit none
        real(kind=dp) :: cluster_count_av, max_size_av, mean_size_av, fract_size_av, av_dist_av
        integer, allocatable, dimension (:) :: max_size, max_size_id
        real(kind=dp), allocatable, dimension (:) :: mean_size, fract_size
	integer :: i_frame, i_cluster, i_mol

	allocate(max_size(end_frame),max_size_id(end_frame),mean_size(end_frame),fract_size(end_frame))

        open(unit = 1, file = 'cluster_out.dat', status = 'unknown')
        100 format (i5, i5, f6.2, i5, f6.2, f6.2)
	200 format (i5, i5, f7.2, 1000i5)
        300 format (a6, f6.2, f6.2, f6.2, f6.2, f6.2)

        ! determine size of largest cluster for all frames
        max_size = maxval(cluster_size,dim=2)
        
        ! determine index of largest cluster for all frames
        max_size_id = maxloc(cluster_size,dim=2)

	write(1,*) "*** all out_style ***"
        write(1,*) "Row 1 Column Numbers for Each Frame:"
	write(1,*) "1. Frame Number"
	write(1,*) "2. Number of Clusters"
	write(1,*) "3. Mean Cluster Size"
	write(1,*) "4. Maximum Cluster Size"
	write(1,*) "5. Fraction of Molecules in Largest Cluster"
	write(1,*) "6. Average Distance Between Molecules"
	write(1,*) "   1    2     3    4     5     6"
	if (out_style.eq.1) then
		write(1,*)
		write(1,*) "*** out_style = 1 ***"
		write(1,*) "As above, followed by n rows for each frame"
		write(1,*) "n = number of clusters in the frame"
		write(1,*) "Each row has the format i, s, t, mol"
		write(1,*) "i = cluster index"
		write(1,*) "s = number of molecules in cluster"
		write(1,*) "t = existence time of exact cluster"
		write(1,*) "mol = indices of molecules in the cluster"
	end if
	write(1,*) "--------------------------------------------------------------"	
	open(unit=2,file='data.dat',status='unknown')
        ! loop over frames
        do i_frame = start_frame, end_frame, space_frame

                mean_size(i_frame) = real(n_mol) / real(cluster_count(i_frame))
                fract_size(i_frame) = real(max_size(i_frame)) / real(n_mol)
                write(1,100) i_frame,cluster_count(i_frame),mean_size(i_frame),max_size(i_frame),&
			& fract_size(i_frame),av_dist(i_frame)
                write(2,*) real(i_frame)*dt,cluster_count(i_frame),mean_size(i_frame),max_size(i_frame)
		! --- write all clusters
		if (out_style.eq.1) then

                	! --- loop over all clusters in this frame
                	do i_cluster = 1, cluster_count(i_frame)
			
				write(1,200) i_cluster, cluster_size(i_frame, i_cluster),real(cluster_time(i_frame,i_cluster))*dt,&
				& (cluster_index(i_frame,i_cluster,i_mol), i_mol=1,cluster_size(i_frame,i_cluster))
				
			end do

                end if

        end do
	write(1,*) "--------------------------------------------------------------"
        cluster_count_av = real(sum(cluster_count)) / real(n_frame)
        mean_size_av = sum(mean_size) / real(n_frame)
        max_size_av = real(sum(max_size)) / real(n_frame)
        fract_size_av = sum(fract_size) / real(n_frame)
	av_dist_av = sum(av_dist) / real(n_frame)
        write(1,300) "Aver: ", cluster_count_av, mean_size_av, max_size_av, fract_size_av, av_dist_av
	close(unit = 1)

	! --- calculate and write out cluster size distribution	
	open(unit = 2, file='dist_out.dat', status='unknown')
	cluster_size_dist(:) = 0.0d0
	do i_frame = start_frame, end_frame, space_frame
		do i_cluster = 1, cluster_count(i_frame)
			do i_mol = 1, n_mol
				if (cluster_size(i_frame,i_cluster).eq.i_mol) then
					cluster_size_dist(i_mol) = cluster_size_dist(i_mol) + 1.0d0
				end if
			end do
		end do		
	end do
	do i_mol = 1, n_mol
		cluster_size_dist(i_mol) = cluster_size_dist(i_mol) / real(n_frame)
		write(2,*) i_mol, cluster_size_dist(i_mol)
	end do
        close(unit = 2)

end subroutine write_out

recursive subroutine mol_cluster(n, m, p)
        use global, only: neighbour, cluster_assign, neigh_num, n_cluster, cluster_index
	implicit none
        integer :: n, m, p, o, i_neigh, n_neigh

        cluster_assign(n, m) = .true.
	cluster_index(n,n_cluster,p) = m 
        n_neigh = neigh_num(n, m)

        if (n_neigh.gt.0) then

                do i_neigh = 1, n_neigh

                        o = neighbour(n, m, i_neigh)
                        if (cluster_assign(n,o)) cycle
                        p = p + 1
                        call mol_cluster(n, o, p)

                end do
        
        end if

end subroutine
