# Bug Report: Index Mismatch in `cal_surface` Function

## Summary

The `cal_surface` function in the BFS calibration code has an index mismatch bug. The function's comment and internal logic expect parameters `[Wb, ALPHA, Ks, Kz]`, but the function is actually called with `[Lb, Wb, ALPHA, Ks]`. This causes incorrect parameter validation checks and may lead to premature termination of the optimization.

## Location

**File:** `source/Rfunctions.bfs_calibration_sub.R`  
**Function:** `cal_surface` (lines 147-162)  
**Calling code:** `source/Rfunction.bfs_calibrate.R` (lines 61-66, 133-136, 169)

## Expected vs Actual Behavior

### Expected (according to function comment):
```r
#X=c(Wb,ALPHA,Ks,Kz)
```

### Actual (what the function is called with):
```r
X=c(Lb,Wb,ALPHA,Ks)
LOGX=log(X,10)
tmp=optim(LOGX,cal_surface,...)
```

## The Bug

### 1. Parameter Validation Check (`bad_X`)

**Location:** Line 150 in `Rfunctions.bfs_calibration_sub.R`

```r
bad_X=c(10^x[1]>basin_char[1]/basin_char[2],x[2]>(-1),x[3]<(-8),x[3]>5,x[4]<(-8),x[4]>5)
```

**Problem:** The function expects `[Wb,ALPHA,Ks,Kz]` but receives `[Lb,Wb,ALPHA,Ks]`, causing:
- `x[1]` expects `Wb` (log10) but gets `Lb` (log10) 
- `x[2]` expects `ALPHA` (log10) but gets `Wb` (log10) ← **This is the critical bug**
- `x[3]` expects `Ks` (log10) but gets `ALPHA` (log10)
- `x[4]` expects `Kz` (log10) but gets `Ks` (log10)

**Impact:** The check `x[2] > -1` is intended to validate that `ALPHA > 0.1` (since log10(0.1) = -1), but instead it checks `Wb > -1`. Since `Wb` (log10) is typically around 3.29, this check is always `TRUE`, causing the function to return the penalty value of 100.0.

### 2. Function Body Parameter Usage

**Location:** Lines 153-157 in `Rfunctions.bfs_calibration_sub.R`

```r
wb=10^x[1]
basin_char[4]=wb

a=10^x[2];ks=10^x[3];kz=10^x[4]
gw_hyd[1]=a;gw_hyd[3]=ks;gw_hyd[5]=kz
```

**Problem:** The function body uses the same wrong indices:
- `wb=10^x[1]` expects `Wb` but gets `Lb`
- `a=10^x[2]` expects `ALPHA` but gets `Wb`
- `ks=10^x[3]` expects `Ks` but gets `ALPHA`
- `kz=10^x[4]` expects `Kz` but gets `Ks`

**Note:** This bug in the function body is less critical because R uses `tmp$par` from `optim()` (line 69 in `Rfunction.bfs_calibrate.R`), not the values calculated inside `cal_surface`. However, it does affect the objective function calculation, which may cause optimization issues.

## Example of the Bug in Action

### Test Case:
```r
# After cal_initial optimization:
Lb = 10590.707791
Wb = 1951.818548
ALPHA = 0.009451
Ks = 2309.369646

X = c(Lb, Wb, ALPHA, Ks)
LOGX = log(X, 10)
# LOGX = [4.024925, 3.290439, -2.024526, 3.363493]
```

### What R's `cal_surface` does:
```r
# R expects [Wb,ALPHA,Ks,Kz] but gets [Lb,Wb,ALPHA,Ks]
# So in cal_surface:
x[1] = 4.024925 (log10(Lb)) - R expects log10(Wb) ✗
x[2] = 3.290439 (log10(Wb)) - R expects log10(ALPHA) ✗
x[3] = -2.024526 (log10(ALPHA)) - R expects log10(Ks) ✗
x[4] = 3.363493 (log10(Ks)) - R expects log10(Kz) ✗

# R's check: x[2] > -1
# R expects: ALPHA (log10) > -1, i.e., ALPHA > 0.1
# But R actually checks: Wb (log10) > -1
# Result: 3.290439 > -1 = TRUE → returns 100.0
```

## Impact

1. **Optimization Performance:** The function returns 100.0 (penalty value) for most parameter combinations, which may cause the optimization to:
   - Terminate early (after only 5 iterations in observed cases)
   - Return the initial guess instead of optimized parameters
   - Fail to properly optimize surface parameters

2. **Calibration Results:** While R still produces valid calibration results (because it uses `tmp$par` from `optim()`), the optimization may be less effective than intended.

3. **Code Correctness:** The function does not behave as documented, which could lead to confusion and maintenance issues.

## Suggested Fix

### Option 1: Update `cal_surface` to match the actual call signature

Update the function to expect `[Lb, Wb, ALPHA, Ks]`:

```r
cal_surface=function(x,tmp.q,dys,timestep,error_basis,basin_char,gw_hyd,flow){
#CALIBRATE SURFACE
#X=c(Lb,Wb,ALPHA,Ks)  # Updated comment
bad_X=c(10^x[1]>basin_char[1]/basin_char[2],  # 10^Lb > area/lb (skip this check or adjust)
        x[2]>(-1),  # Wb > -1 (or remove this check if not needed)
        x[3]<(-8),  # ALPHA < -8
        x[3]>5,     # ALPHA > 5
        x[4]<(-8),  # Ks < -8
        x[4]>5)     # Ks > 5

if(any(bad_X)){obj=100} else {
  lb=10^x[1]  # Extract Lb (x[1] is Lb in log10)
  wb=10^x[2]  # Extract Wb (x[2] is Wb in log10)
  basin_char[2]=lb  # Update Lb (basin_char[2] is Lb)
  basin_char[4]=wb  # Update Wb (basin_char[4] is Wb)

  a=10^x[3];ks=10^x[4]  # Extract ALPHA and Ks (x[3] is ALPHA, x[4] is Ks in log10)
  gw_hyd[1]=a;gw_hyd[3]=ks  # Update ALPHA (gw_hyd[1]) and Ks (gw_hyd[3])
  # Note: Kz is not updated in cal_surface (gw_hyd[5] remains unchanged)

  out=bfs(tmp.q,dys,timestep,error_basis,basin_char,gw_hyd,flow)
  obj=objective(bfs_out,prec=flow[5])}
  
obj}
```

### Option 2: Update the calling code to match the function signature

Change the calling code to pass `[Wb, ALPHA, Ks, Kz]` instead of `[Lb, Wb, ALPHA, Ks]`:

```r
# In Rfunction.bfs_calibrate.R, change:
X=c(Lb,Wb,ALPHA,Ks)  # Old
# To:
X=c(Wb,ALPHA,Ks,Kz)  # New

LOGX=log(X,10)
tmp=optim(LOGX,cal_surface,...)

# And update parameter extraction:
if(all(is.finite(tmp$par))){
  Wb=10^tmp$par[1];ALPHA=10^tmp$par[2];Ks=10^tmp$par[3];Kz=10^tmp$par[4]
  # Lb is not updated in cal_surface
  ...
}
```

## Recommendation

**Option 1 is recommended** because:
1. The current calling pattern `[Lb, Wb, ALPHA, Ks]` makes sense - it allows `cal_surface` to update `Lb` along with the other parameters
2. The function body logic suggests it was intended to work with `[Lb, Wb, ALPHA, Ks]` (it updates `basin_char[2]` which is `Lb`)
3. Only the comment and `bad_X` check need updating

## Verification

To verify the bug, run the calibration with diagnostics and observe:
- `cal_surface` returns 100.0 after only 5 iterations
- The optimization terminates early
- Parameters may not be optimally adjusted

After fixing, `cal_surface` should:
- Return proper objective function values (not 100.0)
- Run for more iterations
- Produce better optimized parameters

## Additional Notes

- The bug does not prevent R from producing valid calibration results because R uses `tmp$par` from `optim()`, not the values calculated inside `cal_surface`
- However, the optimization is likely less effective than intended
- The Python implementation has been updated to match R's buggy behavior for consistency, but this should be fixed in the R code

## Contact

If you have questions or need clarification, please refer to the diagnostic output from `run_calibration_12167000_diagnostics.R` which shows the bug in action.

