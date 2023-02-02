module global

        use iso_fortran_env, only : int64
        integer, parameter :: dp = SELECTED_REAL_KIND(15,99)

        ! --- global variables
        integer :: n_frame, start_frame, end_frame, space_frame, n_atm, n_mol, n_mol_atm
        integer :: i_frame, i_trj, i_typ, i_atm, i_mol, j_mol, i_mol_atm
        real(kind=dp) :: mol_mas, dt
        real(kind=dp), allocatable, dimension (:) :: av_dist, atm_mas, mass
        real(kind=dp), allocatable, dimension (:,:,:) :: x, y, z, dist
        real(kind=dp), allocatable, dimension (:) :: lx, ly, lz
        character(len=2), allocatable, dimension (:) :: type, atm_typ
        real(kind=dp), allocatable, dimension (:,:) :: x_mol, y_mol, z_mol
        real(kind=dp) :: theta_x, theta_y, theta_z, theta_x_sum, theta_y_sum, theta_z_sum
	real(kind=dp) :: xi_x, xi_y, xi_z, xi_x_sum, xi_y_sum, xi_z_sum
        real(kind=dp) :: zeta_x, zeta_y, zeta_z, zeta_x_sum, zeta_y_sum, zeta_z_sum
	real(kind=dp), parameter :: pi = 3.14159265359d0 
	real(kind=dp) :: rdf_max
	integer :: n_bin

	! --- cluster analysis variables
	integer :: neigh_count, n_cluster, size_count, out_style
        integer, allocatable, dimension (:,:,:) :: neighbour, cluster_index
        integer, allocatable, dimension (:,:) :: neigh_num, cluster_size, cluster_time
        integer, allocatable, dimension (:) :: cluster_count
	real(kind=dp), allocatable, dimension (:) :: cluster_size_dist
        logical, allocatable, dimension (:,:) :: cluster_assign
        logical :: all_cluster
	real(kind=dp) :: cluster_cut
	integer, parameter :: max_neighbour = 40

        namelist /input/ start_frame,end_frame,space_frame,dt,n_mol,n_typ,n_bin,rdf_max,&
			& cluster_cut,out_style
	namelist /atoms/ type,mass

end module global

