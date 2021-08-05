function [EKF_track]=EKF_3stateUnknownData(EKF_track,I_P,Q_P)
% Performs carrier tracking via an EKF without data wipeoff
%
%[EKF_track] = EKF_3stateUnknownData(EKF_track,I_P,Q_P)
%
%   Inputs:
%       EKF_track       - EKF parameters (state estimates, covariance
%       matrices)
%       I_P             - In-phase accumulation
%       Q_P             - Quadrature accumulation
%   Outputs:
%       EKF_track       -EKF parameters (state estimates, covariance
%       matrices)
%
% Written by P Blunt 2015
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------
% Propagate the state vector estimate
EKF_track.X_minus=EKF_track.Phi*EKF_track.X_plus;
% Propagate the error covariance matrix 
EKF_track.P_minus=EKF_track.Phi*EKF_track.P_plus*(EKF_track.Phi')+EKF_track.Q;
% Calculate the Kalman gain matrix
EKF_track.K=EKF_track.P_minus*(EKF_track.H')/(EKF_track.H*EKF_track.P_minus*(EKF_track.H')+EKF_track.R);
% Implement carrier loop discriminator (phase detector)
Patan = atan(Q_P/I_P);
% Carrier phase discriminator normalisation function
Natan = 1;
% Formulate the measurement
Z_k = Natan * Patan;
% Update the state vector estimate 
EKF_track.X_plus= EKF_track.X_minus+EKF_track.K*(Z_k -(EKF_track.H*EKF_track.X_minus));
% Update the error covariance matrix 
EKF_track.P_plus=(EKF_track.I-EKF_track.K*EKF_track.H)*EKF_track.P_minus;
                   
end

