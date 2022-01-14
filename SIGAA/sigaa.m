warning('off');
clear; clc;
dbstop if error
disp('Method 1: Use local .bmp images and .INF file to create input excel and study it.')
disp('Method 2: Use previously created input excel and study it.')
disp('Method 3: Use local .bmp images and .INF file to create input excel.')
mode = input('What Method you wish to use 1, 2 or 3? ')
clearAllMemoizedCaches


addpath('./lib') 

switch mode
    case 1
    %main program for Calcium
    
    localdir = pwd
    [sheet testname] = patternfindercircle(localdir);
    oscillations(sheet,testname)
    
    case 2

    localdir = pwd;
    testname=input('Name the test you want to study.','s');    
    oscillations([localdir '\outputs\' testname '\ca2Function.xlsx'],[localdir '\outputs\' testname])
    
    case 3

    localdir = pwd;
    
    [excelname testname] = patternfindercircle(localdir)
     
end
close all;
disp('Execution finished.')   
