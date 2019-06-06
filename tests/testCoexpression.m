function test_suite=testCoexpression()
        initTestSuite;
end

function testCoexpressionSimple()
	pkg load nan;
        assertTrue(isequal(Coexpression([1 1 1;1 1 1]),[1 0 0;0 1 0;0 0 1]))
end
