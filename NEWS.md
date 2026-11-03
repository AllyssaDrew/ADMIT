# ADMIT 0.3.2

* A plateau that reaches an end of the search grid is now marked as censored there (`plateau_censored`, and
  `<=` / `+` when printed). The grid's lower bound is `max(min(q_hat > 0), 1e-5)`, so on a cohort whose
  participants reach very low global ancestry the lower end of the plateau is usually the grid's edge, and
  printing it as a bare number invites it into a paper as though someone had measured it.

