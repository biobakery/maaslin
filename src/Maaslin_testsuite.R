#!/usr/bin/env Rscript

# This is a driver script to execute the suite of Maaslin tests.
# To run: ./Maaslin_testsuite.R

library('testthat')

# Find all of the test folders
working_directory=getwd()
test_folders<-dir(path=working_directory,pattern="test-")

for (test_folder in test_folders) 
{
    # Find the name of the R file to source
    test=gsub("test-","",test_folder)
    print(paste("Running tests for Maaslin:",test))
    R_source=file.path(working_directory,"lib",paste(test,".R",sep=""))
    
    # Source the file to test
    if (file.exists(R_source))
    {
        source(R_source,chdir=T)
    }

    # Run all tests in the folder
    test_dir(test_folder,reporter='Summary')
}

