{
  inputs,
  cell,
}: let
  inherit (cell.packages) epsteinlib epsteinlib_opt;
in {
  inherit epsteinlib epsteinlib_opt;
  default = epsteinlib_opt;
}
