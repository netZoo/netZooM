function test_suite=testCoexpression()
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;
end

function testCoexpressionSimple()
    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    if isOctave
        pkg load nan;
    end
    assertTrue(isequal(Coexpression([1 1 1;1 1 1]),[1 0 0;0 1 0;0 0 1]))
end
