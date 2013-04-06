/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  BFGS tools Copyright (C) 2013  Wilmer Henao
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

bool themin(double * x, double * gradpoints, int i, const int n){
  // This function calculates the minimum coordinate distance from x to gradpoints
  // in a certain gradpoints segment
  double valmin = 1.0;
  for (int j; j < n; j++){
    valmin MIN(valmin, std::abs(x[j] - gradpoints[j + (i * n)]));
  }
  if (valmin < 1e-14)
    return false
  else
    return true;
}
