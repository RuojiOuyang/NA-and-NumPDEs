object = EquationSolver.cpp EquationSolver.h
SHELL = cmd.exe
B = problem_B.cpp
C = problem_C.cpp
D = problem_D.cpp
E = problem_E.cpp
F = problem_F.cpp

problem_B: $(B) $(object)
	g++ -o $@ $(B)

problem_C: $(C) $(object)
	g++ -o $@ $(C)
	
problem_D: $(D) $(object)
	g++ -o $@ $(D)
	
problem_E: $(E) $(object)
	g++ -o $@ $(E)

problem_F: $(F) $(object)
	g++ -o $@ $(F)

clean:
	del *.exe