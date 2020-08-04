#!/bin/bash

export event_code="CMTSOLUTION_201305240544A_GCMT"

#python step0_REAL.py CMTSOLUTION_201305240544A_GCMT

read -p "Proceed to step 1 (y/n)? " answer
case ${answer:0:1} in
    y|Y )
       python step1_REAL_check_data.py CMTSOLUTION_201305240544A_GCMT
       python step2_REAL_flip.py CMTSOLUTION_201305240544A_GCMT
       python automatic_time_windowing_FILTERING_TRIAL.py CMTSOLUTION_201305240544A_GCMT
       python ../derivatives_calculation/calculate_derivatives.py CMTSOLUTION_201305240544A_GCMT
       python cut_and_filter_traces.py CMTSOLUTION_201305240544A_GCMT

       read -p "Proceed to final check (y/n)? " answer2
       case ${answer:0:1} in 
       y|Y )
        python final_check.py CMTSOLUTION_201305240544A_GCMT
       ;;
       * )
        exit 
        ;;
       esac  

    ;;
    * )
        exit
    ;;
esac
