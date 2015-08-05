% Plots the computed fields
close all
clear all
clc

path = '../output/';

% Load data
[E, Ne] = zmread(strcat(path, 'E.txt'));
[position, Np] = zmread(strcat(path, 'position.txt'));
[shape, Ns] = imread(strcat(path, 'shape.txt'));

clear path

% Prepare data
N = Np(1);
M = Np(2);

if (N ~= Ns(1) && M ~= Ns(2) && N*M ~= Ne)
    error('Dimensions of matrices do not match!');
end

Et = reshape(E, N, M);
Ea = 20.0*log10(abs(Et / max(Et(:))));

xlimits = [round(min(real(position(:)))) round(max(real(position(:))))];
ylimits = [round(min(imag(position(:)))) round(max(imag(position(:))))];

X = reshape(real(position), N, M);
Y = reshape(imag(position), N, M);

plot_shape(shape, X, Y, xlimits, ylimits, Ea);