clc;clear all;close all
tic
run('main_area1.m');
run('main_area2.m');
totalTime = toc;
disp(['Total running time: ' num2str(totalTime) ' seconds']);