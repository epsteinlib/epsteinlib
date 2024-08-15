# SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
#
# SPDX-License-Identifier: AGPL-3.0-only
{
  inputs,
  cell,
}: let
  inherit (cell.packages) epsteinlib epsteinlib_dbg;
in {
  inherit epsteinlib epsteinlib_dbg;
  default = epsteinlib;
}
