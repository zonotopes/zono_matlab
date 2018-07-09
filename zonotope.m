function volume = zonotope(star, computeAdditionalStatistics)
% This function computes the zonotope from a given set of generators,
% stored row-by-row in the input matrix 'star'.
%
% When the input argument computeAdditionalStatistics is passed to the
% function and is true, the function computes additional statistics
% and display them on the screen.
% Such statistics are the one useful in economics, ad described
% in the following paper:
%
% Dosi G, Grazzi M, Marengo L, and Settepanella S,
% "Production Theory: Accounting for Firm Heterogeneity and Technical Change",
% The Journal of Industrial Economics}, Vol. 64, No. 4, pp. 875--907, 2016,
% doi:10.1111/joie.12128, http://dx.doi.org/10.1111/joie.12128
%
% % Example
%
% fprintf('Given the following matrix of ten 3D generators\n\n');
%
% star = [ 9.7974838e-01   6.2406009e-01   6.0986665e-01 ;
%          4.3886997e-01   6.7913554e-01   6.1766639e-01 ;
%          1.1111922e-01   3.9551522e-01   8.5944231e-01 ;
%          2.5806470e-01   3.6743665e-01   8.0548942e-01 ;
%          4.0871985e-01   9.8798200e-01   5.7672152e-01 ;
%          5.9489607e-01   3.7738866e-02   1.8292247e-01 ;
%          2.6221175e-01   8.8516801e-01   2.3993201e-01 ;
%          6.0284309e-01   9.1328683e-01   8.8651193e-01 ;
%          7.1121578e-01   7.9618387e-01   2.8674152e-02 ;
%          2.2174673e-01   9.8712279e-02   4.8990139e-01 ];
%
% disp(star);
% vol = zonotope(star);
% fprintf('The volume of the associated zonotope is: %g\n', vol);
% fprintf('(the known correct volume for this dataset is: 20.1948)\n');
%

%-----------------------------------------------------------------------
%   ZONOTOPE LIBRARY ver 1.2
%
%   Copyright (c) 2017, Marco Cococcioni.
%   All rights reserved.
%
%   Redistribution and use in source and binary forms, with or without
%   modification, are permitted provided that the following conditions are met:
%
%   1. Redistributions of source code must retain the above copyright notice, this
%      list of conditions and the following disclaimer.
%   2. Redistributions in binary form must reproduce the above copyright notice,
%      this list of conditions and the following disclaimer in the documentation
%      and/or other materials provided with the distribution.
%
%   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%   ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
%   ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%-----------------------------------------------------------------------
if nargin < 2
    computeAdditionalStatistics = true;
end

if computeAdditionalStatistics
    fprintf('-----------------------------------------------------------------------\n');
    fprintf('ZONOTOPE LIBRARY VER 1.2                       \n');
    fprintf('-----------------------------------------------------------------------\n');
    fprintf('INPUT: SET OF GENERATORS             \n');
    fprintf('N. of dimensions (including the output): %d\n',size(star,2));
    fprintf('N. of generators: %d\n', size(star,1));
    fprintf('-----------------------------------------------------------------------\n');
    fprintf('... the computation of the volume has started (it can take a while) ...\n');
end
  

tic % starts the times   % bt = cputime; % begin time for building the zonotope

dim = size(star, 2); % dimensions of the zonotope
n = size(star, 1);  % number of generators

if ( n < dim )
    error('Zonotope::Error: the number of generators must be greater than the dimension of the space');
end

face = zeros(n, 1); % vectors lying on the solid angle at the origin of the zonohedron

counters = zeros(dim,1);
for i = 1:dim-1
    counters(i) = i;
end
counters(dim) = n + 1;

a = zeros(dim, dim);

volume = 0; % let us initialize the volume to be computed to zero

while_counter = 0; % number of while iteration of the next while)
keep_going = true; % auxiliary boolean variable that control the exit from the while
issue_linear_dependency_warning = true;

while (keep_going) % MAIN LOOP NEEDED TO COMPUTE THE ZONOTOPE    
    a(1:dim-1, 1:dim) = star(counters(1:dim-1), 1:dim);    
    a(dim,1:dim) = star(counters(1),1:dim);
    
    % Compute the normal of the vectors in a
    % A normal can be obtained by computing the null space of a
    % (i.e., its kernel)
        
    normal = null(a);    % normal = kernel(a);
    
    if (size(normal,2) > 1)
        if ( issue_linear_dependency_warning )
            warning('Some of the generators are linearly dependent. The volume might be inaccurate')
            issue_linear_dependency_warning = false;
        end
        normal = normal(:,1); %take first column only
    end
    
    normal = normal/norm(normal);  % normal = normalize(normal);
    
    % I need to find orientation of the normal
    a(dim, 1:dim) = normal(1:dim);
    s = det(a);
    
    if s < 0
        normal = -normal;
    end
    
    aux = sign(star*normal);
    
    starting = 1;
    for i=1:dim
        ending = counters(i);
        K = starting:ending-1;
        face(K) = aux(K);
        
        if ending < n
            face(ending) = 0;
        end
        starting = counters(i) + 1;
    end
    
    % base = project(face, star, dim);
    base = sum(star.*face(:,ones(1,dim)));    
    
    % next row computes the scalar product between vector 'base'
    % (a row vector) and normal*abs(s), a column vector
    pyramid = base*normal*abs(s);
    
    % next row updates the volume
    volume = volume + pyramid;
    
    % next row determines if we are done or not
    [ keep_going, counters ] = increase(counters, dim-2, dim-1, size(star,1));
    
    % next row update the number of iterations performed
    while_counter = while_counter + 1;
    
end % MAIN LOOP

% final correction to the volume
volume = volume / (factorial(dim)/2); %they are all N dimensional cones. (2, because we are computing half the faces, now)

%et = cputime; % end time for building the zonotope

% Now compute, if requested, additional statistics other than the volume
if computeAdditionalStatistics,    
%    % Computing tagi, which is the tangent between each input generator
%    %   and the (dim-1)-dimensional input space )
%    tagi = zeros(n,1);    
%    for g = 1:n
%        tagi(g) = star(g,end)/sqrt(sum(star(g,1:end-1).^2));
%    end
    
    fprintf('-----------------------------------------------------------------------\n');
    fprintf('OUTPUT STATISTICS                              \n');
    fprintf('S1: Total volume:  %g\n', volume);    
    
    diagonal = sum(star);
    fprintf('S2: Diagonal: %s\n', num2str(diagonal));
    
    diagonal_norm = sqrt(sum(diagonal.^2));
    fprintf('S3: Diagonal norm: %g\n', diagonal_norm);
    
    total_length_squared = sum(star(:).^2);
    fprintf('S4: Sum of squared norms: %g\n', total_length_squared);
    
    % Compute the Gini index (it depends on volume and diagonal)
    gini = volume / prod(diagonal);
    fprintf('S5: Gini index: %g\n', gini);
    
    tang_diag_input = diagonal(end)/sqrt(sum(diagonal(1:end-1).^2));
    fprintf('S6: Tangent of angle btw. diagonal and the input plane: %g\n', tang_diag_input);
    
    cos1 = diagonal(end) / norm(diagonal);
    fprintf('S7: Cosine against output: %g\n', cos1 );
    
    cos2 = diagonal(1)/sqrt(sum(diagonal(1:end-1).^2));
    fprintf('S8: Cosine of projection of diagonal on input plane with x axis: %g\n', cos2);    
    
    fprintf('S9: Volume against the cube of the norm of the diagonal: %g\n' , volume/(diagonal_norm^3) );
                    
    % ADDITIONAL_STATISTICS (for 3D problems only)
%     if dim == 3
%         fprintf('-----------------------------------------------------------------------\n');
%         fprintf('ADDITIONAL STATISTICS (for 3D problems only)           \n');
        
        % % Compute the Gini index (it depends on volume and diagonal)
        % gini = volume / (diagonal(1) * diagonal(2) * diagonal(3));
        % fprintf('B1: Gini 3D index: %g\n', gini);
        
        %tang_diag_input = diagonal(3)/sqrt(diagonal(1)^2 + diagonal(2)^2);
        %fprintf('B2: Tangent of angle btw. diagonal and the x-y plane: %g\n', tang_diag_input);
        
        %cos1 = diagonal(3) / norm(diagonal);
        %fprintf('B3: Cosine against output: %g\n', cos1 );
        
        %cos2 = diagonal(1)/sqrt(diagonal(1)*diagonal(1) + diagonal(2)*diagonal(2));
        %fprintf('B4: Cosine of projection of diagonal on plane x-y with x axis: %g\n', cos2);    
    
        %fprintf('B5: Volume against the cube of the norm of the diagonal: %g\n' , volume/(diagonal_norm^3) );
    
%        aux = (total_length_squared/3.0).^1.5;
%        mistery_number = binomial(size(star,1),3) * aux;
%        fprintf('B6: Mystery number: %g\n', mistery_number);       
        
%        fprintf('B7: Tangent against input axes: %g\n', diagonal(2)/diagonal(1) );               
%        fprintf('B8: Tangent against input axes: %g\n', xx);        
%        fprintf('B9: Solid angle: %g\n', solidAngle(star, edges));        
%        fprintf('B10: Normalized vectors volume: %g\n', xxx);                
%        fprintf('B11: Volume against diagonal cubed of boundary vectors: %g \n', xxxx);        
        
%     end % 3D specific statistics
    fprintf('-----------------------------------------------------------------------\n');
    %fprintf('Elapsed Time\n%g min\n', (et-bt)/60);
    %fprintf('Elapsed Time\n');
    toc
    fprintf('-----------------------------------------------------------------------\n');
end % END of additional statistics
end % END of zonotope main function

    

%------------------------------------------------------------------
% UTILITY FUNCTIONS
%------------------------------------------------------------------

function [res, counters ] = increase(counters, pos, dim, n)
if ( pos < 0 )
    res = false;
    return
end
if ( counters(pos+1) == ( n - (dim - pos) + 1 ) ),
    [res, counters ] = increase(counters, pos-1, dim, n);
    return
end
counters(pos+1) = counters(pos+1) + 1;
for i=pos+1:dim-1
    counters(i+1) = counters(pos+1)+(i-pos);
end
res = true;
return
end