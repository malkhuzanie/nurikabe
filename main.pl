% ----------------- BEGIN UTILS ----------------- %
%
colour_code(red, "\e[31m").
colour_code(green, "\e[32m").
colour_code(blue, "\e[34m").
colour_code(white, "\e[37m").
colour_code(reset, "\e[0m").


colour(b, _colour) :- colour_code(blue, _colour), !.
colour('E', _colour) :- colour_code(white, _colour), !.
colour(_, _colour) :- colour_code(green, _colour), !.

print_grid(_, _j, _, M) :- _j >= M, write("|"), !.
print_grid(_i, _j, N, M) :- 
  fxd_cell(_i, _j, C), 
  colour(C, _colour),
  colour_code(reset, _r_colour),
  format("| ~w~w~w ", [_colour, C, _r_colour]), 
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
  at(_grid, _i, _j, C), 
  colour(C, _colour),
  colour_code(reset, _r_colour),
  format("| ~w~w~w ", [_colour, C, _r_colour]), 
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


% ------------------ END UTILS ------------------ %

:- include('utils.pl').

/*
 * fxd_cell represents a cell located at (Row, Col)
 * which have a number of islands = Num
 */

:- dynamic fxd_cell/3.

/*
 * solve_cell
 * colour a cell located at (Row, Col)
 * with the colour = Num
 */

solve_cell(Row, Col, Num) :-
  asserta(fxd_cell(Row, Col, Num)).


% contains the list grid representation
:- dynamic grid/1.

/* each cell in the grid is either a
 * Num if it contains a number
 * green | g for a green cell
 * blue | b for a blue cell
 * E for empty cell
 */

:- dynamic row_count/1.
:- dynamic col_count/1.

read_dimensions() :-
  write("Enter the number of rows: "),
  read(N),
  asserta(row_count(N)),
  write("Enter the number of columns: "),
  read(M),
  asserta(col_count(M)).

init_grid() :- 
  read_dimensions(),
  row_count(N),
  col_count(M),
  create_grid(0, 0, N, M, 'E').

main :- 
  clear(),
  init_grid().

clear :-
  retractall(fxd_cell(_, _, _)),
  retractall(row_count(_)),
  retractall(col_count(_)).

% check if a cell is outside the borders of the grid
outside(_i, _j) :-
  row_count(N), 
  col_count(M),
  (_i < 0 ; _i >= N ; _j < 0 ; _j >= M).

% check if a cell is inside the borders of the grid
inside(_i, _j) :-
  row_count(N), 
  col_count(M),
  (_i >= 0, _i < N, _j >= 0, _j < M).

% check if a cell is empty 
is_empty(_i, _j, 1) :- fxd_cell(_i, _j, 'E'), !.
is_empty(_, _, 0).

% print grid 
print() :- 
  row_count(N),
  col_count(M),
  print(0, 0, N, M).

:- dynamic vis/2.

clear_vis :- retractall(vis(_, _)).

dfs(_i, _j, _num, 0) :- 
  (\+ fxd_cell(_i, _j, _num) ; vis(_i, _j)), !.
dfs(_i, _j, _, 0) :- outside(_i, _j), !.
dfs(_i, _j, _num, _count) :- 
  asserta(vis(_i, _j)),
  (L is _j - 1, (dfs(_i, L, _num, _lc) -> true ; _lc = 0)),
  (U is _i - 1, (dfs(U, _j, _num, _uc) -> true ; _uc = 0)),
  (R is _j + 1, (dfs(_i, R, _num, _rc) -> true ; _rc = 0)),
  (D is _i + 1, (dfs(D, _j, _num, _dc) -> true ; _dc = 0)),
  _count is 1 + _lc + _uc + _rc + _dc.

exists(_i, _j, _num, _out) :- 
  (fxd_cell(_i, _j, _num) -> _out = 1 ; _out = 0).

one_sea :- 
  once(fxd_cell(X, Y, b)),
  row_count(N),
  count(0, 0, N, b, _s_c_count),
  clear_vis(),
  dfs(X, Y, b, _count),
  _count = _s_c_count.

is_2_by_2_sea(_i, _j, 1) :-
  (
    (fxd_cell(_i, _j, b)),
    (R is _j + 1, fxd_cell(_i, R, b)),
    (D is _i + 1, fxd_cell(D, _j, b)),
    (IS is _i + 1, JS is _j + 1, fxd_cell(IS, JS, b))
  ), !.
is_2_by_2_sea(_, _, 0).

count_row_2_by_2_sea(_i, _j, 0) :- outside(_i, _j), !.
count_row_2_by_2_sea(_i, _j, _count) :-
  is_2_by_2_sea(_i, _j, _a),
  J is _j + 1,
  count_row_2_by_2_sea(_i, J, _b),
  _count is _a + _b.

count_2_by_2_sea(_i, _j, 0) :- outside(_i, _j), !.
count_2_by_2_sea(_i, _j, _count) :- 
  count_row_2_by_2_sea(_i, _j, _a),
  I is _i + 1,
  count_2_by_2_sea(I, _j, _b),
  _count is _a + _b.

no_2_by_2_sea :- 
  count_2_by_2_sea(0, 0, _count), 
  _count = 0.

is_numbered_cell(_i, _j) :-
  fxd_cell(_i, _j, _v),
  _v \= 'E',
  _v \= g, 
  _v \= b.

is_empty_cell(_i, _j) :- fxd_cell(_i, _j, 'E').

traverse_island(_i, _j, 0) :- outside(_i, _j), !.
traverse_island(_i, _j, 0) :- 
  (fxd_cell(_i, _j, b) ; vis(_i, _j)), !.
traverse_island(_i, _j, _size) :- 
  asserta(vis(_i, _j)),
  (L is _j - 1, (traverse_island(_i, L, _lc) -> true ; _lc = 0)),
  (U is _i - 1, (traverse_island(U, _j, _uc) -> true ; _uc = 0)),
  (R is _j + 1, (traverse_island(_i, R, _rc) -> true ; _rc = 0)),
  (D is _i + 1, (traverse_island(D, _j, _dc) -> true ; _dc = 0)),
  _size is 1 + _lc + _uc + _rc + _dc.

island_size_match(_i, _j, _size, _out) :-
  fxd_cell(_i, _j, _num),
  (_num = _size -> _out = 1 ; _out = 0).

island_component_size_h(_i, _j, 0, 0) :- outside(_i, _j), !.
island_component_size_h(_i, _j, A, B) :-
  (is_numbered_cell(_i, _j) 
  ->
    (
      traverse_island(_i, _j, _size),
      island_size_match(_i, _j, _size, _out),
      C is 1
    )
  ;
    (
      C is 0,
      _out is 0
    )
  ),
  J is _j + 1,
  (island_component_size_h(_i, J, NA, NB) -> true ; NA = 0, NB = 0),
  (A is C + NA, B is _out + NB).

island_component_size(_i, _j, 0, 0) :- outside(_i, _j), !.
island_component_size(_i, _j, A, B) :-
  island_component_size_h(_i, _j, RA, RB),
  I is _i + 1, 
  (island_component_size(I, _j, NA, NB) -> true ; NA = 0, NB = 0),
  (A is RA + NA, B is RB + NB).

island_number_equals_size :-
  clear_vis(),
  island_component_size(0, 0, A, B),
  A = B.

% check if each numbered cell belongs to a
% certain island, i.e. a group of connected
% green cells

find_visited_islands(_cells) :- 
  findall((X, Y), (fxd_cell(X, Y, g), vis(X, Y)), _cells).

one_fixed_cell_in_island :- 
  find_visited_islands(_cells),
  row_count(N),
  count(0, 0, N, g, _count),
  length(_cells, _len),
  _count == _len.

is_valid_solution :-
  one_sea, 
  no_2_by_2_sea,
  island_number_equals_size,
  one_fixed_cell_in_island.

% :- include('grid.pl').


% --------------------------------------------------------------- %
% ------------------------> solve section <---------------------- %
% --------------------------------------------------------------- %

set(_i, _j, _) :- outside(_i, _j), !.
set(_i, _j, _num) :-
  once(retract(fxd_cell(_i, _j, _)) ; true),
  solve_cell(_i, _j, _num).

set_if(_i, _j, _, _pre) :- 
  (\+ fxd_cell(_i, _j, _pre) ; outside(_i, _j)), !.
set_if(_i, _j, _num, _) :-
  once(retract(fxd_cell(_i, _j, _)) ; true),
  solve_cell(_i, _j, _num).

create_grid(_i, _j, M, _x) :-
  _j < M,
  solve_cell(_i, _j, _x),
  J is _j + 1,
  create_grid(_i, J, M, _x).

create_grid(_i, _j, N, M, _x) :-
  _i < N,
  (once(create_grid(_i, _j, M, _x)) ; true),
  I is _i + 1, 
  (once(create_grid(I, _j, N, M, _x)) ; true).

% ------------------------------------------------- %
% ------------------ Island of 1 ------------------ %
% ------------------------------------------------- %

solve_one(_i, _j) :-
  (L is _j - 1, once(set_if(_i, L, b, 'E')) ; true),
  (U is _i - 1, once(set_if(U, _j, b, 'E')) ; true),
  (R is _j + 1, once(set_if(_i, R, b, 'E')) ; true),
  (D is _i + 1, once(set_if(D, _j, b, 'E')) ; true).

solve_one_island(_, _j, M) :- _j >= M, !.
solve_one_island(_i, _j, M) :-
  (fxd_cell(_i, _j, 1) -> solve_one(_i, _j) ; true),
  J is _j + 1,
  solve_one_island(_i, J, M).

solve_one_island(_i, _, N, _) :- _i >= N, !.
solve_one_island(_i, _j, N, M) :-
  (once(solve_one_island(_i, _j, M)) ; true),
  I is _i + 1, 
  solve_one_island(I, _j, N, M).

% ------------------------------------------------- %
% ---------------- End Island of 1 ---------------- %
% ------------------------------------------------- %
%
% ------------------------------------------------- %
% --------- Clues separated by one square --------- %
% ------------------------------------------------- %

solve_neighbour_row(_i, _j, _dir) :-
  NJ is _j + _dir * 2,
  (is_numbered_cell(_i, NJ)
  -> 
    (
      J is _j + _dir, set_if(_i, J, b, 'E')
    )
  ;
    (
      true
    )
  ).
    
solve_neighbour_col(_i, _j, _dir) :-
  NI is _i + _dir * 2,
  (is_numbered_cell(NI, _j)
  -> 
    (
      I is _i + _dir, set_if(I, _j, b, 'E')
    )
  ;
    (
      true
    )
  ).

solve_neighbour_island(_, _j, M) :- _j >= M, !.
solve_neighbour_island(_i, _j, M) :-
  (is_numbered_cell(_i, _j) 
  -> 
    (
      solve_neighbour_row(_i, _j, +1),
      solve_neighbour_col(_i, _j, +1)
    )
  ; 
    true
  ),
  J is _j + 1,
  solve_neighbour_island(_i, J, M).

solve_neighbour_island(_i, _, N, _) :- _i >= N, !.
solve_neighbour_island(_i, _j, N, M) :-
  (once(solve_neighbour_island(_i, _j, M)) ; true),
  I is _i + 1, 
  solve_neighbour_island(I, _j, N, M).

% ------------------------------------------------- %
% ----------- Diagonally adjacent clues ----------- %
% ------------------------------------------------- %

solve_diag1(_i, _j) :-
  IU is _i - 1, 
  JU is _j + 1,
  (is_numbered_cell(IU, JU) 
  -> 
    (
      (I0 is _i - 1, set_if(I0, _j, b, 'E')),
      (J0 is _j + 1, set_if(_i, J0, b, 'E'))
    )
  ; 
    true
  ).

solve_diag2(_i, _j) :-
  ID is _i + 1, 
  JD is _j + 1,
  (is_numbered_cell(ID, JD) 
  -> 
    (
      (I0 is _i + 1, set_if(I0, _j, b, 'E')),
      (J0 is _j + 1, set_if(_i, J0, b, 'E'))
    )
  ; 
    true
  ).

solve_diag_island(_, _j, M) :- _j >= M, !.
solve_diag_island(_i, _j, M) :-
  (is_numbered_cell(_i, _j) 
  -> 
    (
      solve_diag1(_i, _j), 
      solve_diag2(_i, _j)
    )
  ; 
    true
  ),
  J is _j + 1,
  solve_diag_island(_i, J, M).

solve_diag_island(_i, _, N, _) :- _i >= N, !.
solve_diag_island(_i, _j, N, M) :-
  (once(solve_diag_island(_i, _j, M)) ; true),
  I is _i + 1, 
  solve_diag_island(I, _j, N, M).


% ------------------------------------------------- %
% --------------- Surrounded Sqare ---------------- %
% ------------------------------------------------- %

is_surrounded(_i, _j) :- 
  (L is _j - 1, (outside(_i, L) ; fxd_cell(_i, L, b))),
  (U is _i - 1, (outside(U, _j) ; fxd_cell(U, _j, b))),
  (R is _j + 1, (outside(_i, R) ; fxd_cell(_i, R, b))),
  (D is _i + 1, (outside(D, _j) ; fxd_cell(D, _j, b))).

solve_sur_square(_, _j, M) :- _j >= M, !.
solve_sur_square(_i, _j, M) :-
  (is_empty_cell(_i, _j), is_surrounded(_i, _j)
  -> 
    (
      % write(_i), write(_j),
      set(_i, _j, b)
    )
  ; 
    true
  ),
  J is _j + 1,
  solve_sur_square(_i, J, M).

solve_sur_square(_i, _, N, _) :- _i >= N, !.
solve_sur_square(_i, _j, N, M) :-
  (once(solve_sur_square(_i, _j, M)) ; true),
  I is _i + 1, 
  solve_sur_square(I, _j, N, M).

% ------------------------------------------------- %
% ------------------ Create grid ------------------ %
% ------------------------------------------------- %

create_list([], _, M, M, _) :- !.
create_list([X | Rest], _row, _col, M, _) :-
    _col < M,
    fxd_cell(_row, _col, X),
    _n_col is _col + 1,
    create_list(Rest, _row, _n_col, M, _).

create_fgrid([], N, _, N, _) :- !.
create_fgrid([R | Rest], _row, _col, N, M) :-
    _row < N,
    create_list(R, _row, _col, M, _),
    _n_row is _row + 1,
    create_fgrid(Rest, _n_row, _col, N, M).

% count the number of cells in a grid which have 
% a number equals to _num
count(_i, _j, _, 0) :- outside(_i, _j), !.
count(_i, _j, _num, _count) :-
  (fxd_cell(_i, _j, _num) -> _out = 1 ; _out = 0),
  J is _j + 1,
  count(_i, J, _num, _rem),
  _count is _out + _rem.

count(_i, _, N, _, 0) :- _i >= N, !.
count(_i, _j, N, _num, _count) :-
  count(_i, _j, _num, _r_count),
  I is _i + 1, 
  count(I, _j, N, _num, _rem_count),
  _count is _r_count + _rem_count.

no_empty_cells() :-
  row_count(N),
  count(0, 0, N, _empty_cells),
  _empty_cells = 0.


empty_cells_set(_empty_cells_set) :-
  findall((X, Y), fxd_cell(X, Y, 'E'), _empty_cells_set).

% ------------------------------------------------- %
% ----------------- Find solution ----------------- %
% ------------------------------------------------- %

:- dynamic solution_not_found/1.
:- dynamic solution/1.

solution_not_found(1).

validate_solution :-
  no_empty_cells(),
  is_valid_solution(), 
  solution_not_found(1),
  row_count(N), 
  col_count(M),
  create_fgrid(_grid, 0, 0, N, M), 
  asserta(solution(_grid)),
  retract(solution_not_found(1)),
  write('-------> '),
  write('A valid solution has been found'),
  write(' <-------\n'),
  print(_grid, 0, 0, N, M),
  !.

permutate([]) :- validate_solution(), !.
permutate([(X, Y) | T]) :-
  set(X, Y, b),
  (no_2_by_2_sea() -> (permutate(T) -> true ; true) ; true),
  set(X, Y, g),
  (no_2_by_2_sea() -> (permutate(T) -> true ; true) ; true).

solve :-
  retractall(solution(_)),
  asserta(solution_not_found(1)),
  row_count(N),
  col_count(M),
  (solve_one_island(0, 0, N, M), print()),
  (solve_neighbour_island(0, 0, N, M), print()),
  (solve_diag_island(0, 0, N, M), print()),
  (solve_sur_square(0, 0, N, M), print()),
  empty_cells_set(_set),
  permutate(_set).

