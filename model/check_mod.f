

      MODULE check


        Use param
        Use param1

!                        variables for check_mass_balance
        DOUBLE PRECISION start_time, report_time
        DOUBLE PRECISION accumulation_g, accumulation_X_g(DIM_N_g)
        DOUBLE PRECISION flux_out_g(DIMENSION_BC), flux_in_g(DIMENSION_BC)
        DOUBLE PRECISION flux_out_X_g(DIMENSION_BC,DIM_N_g), flux_in_X_g(DIMENSION_BC,DIM_N_g)
	

      END MODULE check                                                                             
