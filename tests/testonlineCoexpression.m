function test_suite=testonlineCoexpression()
	initTestSuite;
end

function testonlineCoexpressionSimple()
	% Tell if this is Octave (Unit tests) or Matlab
        isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

        % Load statistics package from Octave
        if isOctave
                %we need the nan package because it has a fast implementation of corrcoeff
                %pkg load statistics
            pkg load nan;
        end

       
        for i=1:100
            % 1. correlation for n samples
            x1 =rand(1000,1000);
            c1 = corr(x1);
            m1 = mean(x1);
            std1 = std(x1);

            % 2. correlation without one sample
            x2  = x1(2:end,:);
            c2  = corr(x2);

            % 3. a faster way to do 2 using 1
            [n,m]=size(x1); % n is number of observations
            s1 = x1(1,:); % this is sample i, (it is a row vector)
            % We compute the covariance
            cc1=c1.*(std1'*std1);
            c3=onlineCoexpression(s1,n,m1,std1,cc1);

            % the matrices are comparable
            assert(max(max(abs(c3-c2))) < 1e-12);

        end



end
