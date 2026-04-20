# How to Apply the Fix

## Quick Summary

**File:** `qre_unitary.py`  
**Lines:** 106-107  
**Status:** ✅ Verified with analytical test case

## The Fix

### Before (INCORRECT):
```python
s1 = timestep * config_unitary.error_scale * c1 / eps_trotter
s2 = timestep * math.sqrt(config_unitary.error_scale * c2 / eps_trotter)
```

### After (CORRECT):
```python
s1 = timestep**2 * config_unitary.error_scale * c1 / (2 * eps_trotter)
s2 = (timestep**1.5) * math.sqrt(config_unitary.error_scale * c2 / eps_trotter)
```

## Changes Explained

### Change 1: First-order formula (s1)
- **Add:** `timestep**2` (square the timestep)
- **Add:** Division by `2` 
- **Why:** Matches Equation 145: r = t² * c1 / (2ε)

### Change 2: Second-order formula (s2)
- **Change:** `timestep` → `(timestep**1.5)`
- **Why:** Matches Equation 152: r = t^(3/2) * sqrt(c2 / ε)

## Step-by-Step Instructions

### 1. Make a backup
```bash
cd qre
cp qre_unitary.py qre_unitary.py.backup
```

### 2. Apply the fix

Option A - Manual edit:
- Open `qre_unitary.py` in your editor
- Go to lines 106-107
- Replace the two lines as shown above
- Save the file

Option B - Use sed (Linux/Mac):
```bash
cd qre
sed -i.bak '106s/timestep \* config_unitary.error_scale \* c1 \/ eps_trotter/timestep**2 * config_unitary.error_scale * c1 \/ (2 * eps_trotter)/' qre_unitary.py
sed -i '107s/timestep \* math.sqrt/\(timestep**1.5\) * math.sqrt/' qre_unitary.py
```

### 3. Verify the fix
```bash
python3.11 test_proposed_fix.py
```

You should see:
```
✅ FIX WORKS! (for both first and second order)
✅ ALL TIME POINTS MATCH - FIX IS CORRECT!
```

### 4. Run existing tests
```bash
python3.11 test_correctness.py
```

Should still pass (these test C1/C2 calculation, not step counts).

### 5. Optional: Add comment
Add this comment above lines 106-107:
```python
# Trotter step counts from Childs et al. arXiv:1912.08854v3
# First-order (Eq 145): r = t² * c1 / (2ε)
# Second-order (Eq 152): r = t^(3/2) * sqrt(c2 / ε)
# where ε = eps_trotter = energy_error * timestep / (2π)
s1 = timestep**2 * config_unitary.error_scale * c1 / (2 * eps_trotter)
s2 = (timestep**1.5) * math.sqrt(config_unitary.error_scale * c2 / eps_trotter)
```

## Complete Code Block

Here's the full section with context (lines 104-108):

```python
eps_trotter = config_unitary.energy_error * timestep / (2 * math.pi)
config_general.log(f"-- allowable fractional Trotter error = {eps_trotter}")
# Trotter step counts from Childs et al. arXiv:1912.08854v3
# First-order (Eq 145): r = t² * c1 / (2ε)
# Second-order (Eq 152): r = t^(3/2) * sqrt(c2 / ε)
s1 = timestep**2 * config_unitary.error_scale * c1 / (2 * eps_trotter)
s2 = (timestep**1.5) * math.sqrt(config_unitary.error_scale * c2 / eps_trotter)
config_general.log(f"-- Trotter step count: s1 = {s1}, s2 = {s2}")
```

## What Changes for Users

### Resource Estimates Will Change!

| Scenario | Impact |
|----------|--------|
| First-order, timestep=1.0 | 2× more steps (was undercounting) |
| First-order, timestep=2.0 | 4× more steps |
| First-order, timestep=0.5 | Same steps (by coincidence) |
| Second-order, timestep=1.0 | Same steps (was accidentally correct) |
| Second-order, timestep=2.0 | 1.4× more steps |
| Second-order, timestep=0.5 | 0.7× steps (fewer) |

### If You've Published Results

If you've published resource estimates using:
- First-order Trotter at any timestep
- Second-order Trotter with timestep ≠ 1.0

You should recalculate with the corrected formulas.

## Rollback (if needed)

If something goes wrong:
```bash
cd qre
cp qre_unitary.py.backup qre_unitary.py
```

## Questions?

If the test shows ❌ instead of ✅, or if you have other issues:
1. Check you're editing the right file (`qre_unitary.py` in the `qre` directory)
2. Make sure you're using Python 3.11
3. Verify the changes are exactly as shown
4. Check the git diff to see what changed

Run:
```bash
git diff qre_unitary.py
```

Should show:
```diff
-    s1 = timestep * config_unitary.error_scale * c1 / eps_trotter
-    s2 = timestep * math.sqrt(config_unitary.error_scale * c2 / eps_trotter)
+    s1 = timestep**2 * config_unitary.error_scale * c1 / (2 * eps_trotter)
+    s2 = (timestep**1.5) * math.sqrt(config_unitary.error_scale * c2 / eps_trotter)
```
