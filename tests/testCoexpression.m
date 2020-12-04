function test_suite=testCoexpression()
        initTestSuite;
end

function testCoexpressionSimple()
    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    if isOctave
        pkg load nan;
    end
    assertTrue(isequal(Coexpression([1 1 1;1 1 1]),[1 0 0;0 1 0;0 0 1]))
end
