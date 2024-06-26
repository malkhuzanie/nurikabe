print_grid(_, _j, _, M) :- _j >= M, write("|"), !.
print_grid(_i, _j, N, M) :- 
  % format('(~w, ~w) ', [ _i, _j ]),
  fxd_cell(_i, _j, C), 
  format("| ~w ", C), 
  (J is _j + 1, print_grid(_i, J, N, M)).

print_line(_i, _n, _) :- _i >= _n, !.
print_line(_i, _n, _char) :-
write(_char),
  I is _i + 1, 
  print_line(I, _n, _char).

print(_i, _, N, M) :- 
  (
    _m is 4 * M, 
    print_line(0, _m, '-'), 
    write('-\n')
  ),
  _i >= N, !.
print(_i, _j, N, M) :- 
  (print_grid(_i, _j, N, M), write('\n')), 
  I is _i + 1, 
  print(I, _j, N, M).

at(_grid, _i, _j, _cell) :- 
  nth0(_i, _grid, _row),
  nth0(_j, _row, _cell).


print_grid(_, _, _j, _, M) :- _j >= M, write("|"), !.
print_grid(_grid, _i, _j, N, M) :- 
  % format('(~w, ~w) ', [ _i, _j ]),
  at(_grid, _i, _j, C), 
  format("| ~w ", C), 
  (J is _j + 1, print_grid(_grid, _i, J, N, M)).

print(_, _i, _, N, M) :- 
  (
    _m is 4 * M, 
    print_line(0, _m, '-'), 
    write('-\n')
  ),
  _i >= N, !.
print(_grid, _i, _j, N, M) :- 
  (print_grid(_grid, _i, _j, N, M), write('\n')), 
  I is _i + 1, 
  print(_grid, I, _j, N, M).
