
build_model_newb:
	tr -d '\r' < input_for_mobility.txt > new_input.txt
	gcc -c 2degmobilityUnix_newb.c
	gcc -c inputfile_read_subroutines.c
	gcc -c 2DEGchi_adaptive_subroutine.c
	gcc -c formfactor2deg_subroutine.c
	gcc -o 2degmobilityUnix_newb 2degmobilityUnix_newb.o 2DEGchi_adaptive_subroutine.o inputfile_read_subroutines.o formfactor2deg_subroutine.o -lm
