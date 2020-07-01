        program read_data
c       This program hopefully reads the parameter from the file

            parameter     ( c = 1 )
            real Mo
            open (unit=1,file='scalar_moment.txt')
            read(1,*) Mo
            open (unit = 1, file = )
            area = (Amax * sqrt(3.))* (Amin * sqrt(3.))
            root_area = sqrt (area)
            stress_drop = (c*Mo) / root_area
            print *, area 
            print *, root_area
            print *, stress_drop


        end
