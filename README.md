# IC_visualisation
Routines to extract and plot data produced by the IC dynamics code
scipy.io.FortranFile can read an output file produced by Fortran. Each sequence of writing is read, line after line. 


=====

Default folder for data is ./OUT/

To run the code:
>> python IC_read_data.py
This will run what is in the __main__ part of the file. Functions can
be called separately by `import IC_read_data` and calling
IC_read_data.function() for a given function(). 
TODO: create a separate file which call the functions and do the
extraction and data analysis. (for now, it is in the __main__)

IC_read_data.py : functions to read file + do basic data modification 

IC_plot.py : functions to render figures from the data. 


