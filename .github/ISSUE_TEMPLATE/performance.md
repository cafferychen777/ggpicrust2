---
name: ⚡ Performance Issue
about: Report slow performance, memory issues, or computational problems
title: '[PERFORMANCE] '
labels: 'performance'
assignees: ''

---

## ⚡ Performance Issue

### 🐌 Type of Performance Problem
- [ ] Function runs very slowly
- [ ] High memory usage/out of memory errors
- [ ] Process hangs/freezes
- [ ] Results take much longer than expected
- [ ] Memory usage keeps growing during execution
- [ ] Other: ___________

### 🔍 Affected Function(s)
**Which function(s) are experiencing performance issues?**
- [ ] `pathway_daa()`
- [ ] `pathway_errorbar()`  
- [ ] `pathway_heatmap()`
- [ ] `pathway_pca()`
- [ ] `ggpicrust2()`
- [ ] Data processing functions
- [ ] Multiple functions
- [ ] Other: ___________

### 📊 Data Characteristics
**Your dataset size:**
- Number of samples: 
- Number of features/pathways: 
- Number of metadata variables: 
- File size of abundance data: 
- Total memory usage: 

### ⏱️ Performance Details
**Current performance:**
- Time to complete: 
- Memory usage: 
- CPU usage: 

**Expected performance:**
- Expected time: 
- Comparison with similar analyses: 

### 💻 System Specifications
- **Operating System:** 
- **RAM:** 
- **CPU:** 
- **R Version:** 
- **ggpicrust2 Version:** 

### 📋 Reproducible Example
Please provide a minimal example that demonstrates the performance issue:

```r
# Minimal example showing the performance problem
# Include system.time() or profiling information if possible
library(ggpicrust2)

system.time({
  # Your slow code here
})
```

### 🔧 Analysis Parameters
**What parameters are you using?**
```r
# Please share the parameters you're using that might affect performance
# For example:
method = "..."
p_adjust = "..."
# etc.
```

### 📈 Performance Comparison
**Performance across different conditions:**
- [ ] Same dataset, different parameters
- [ ] Different dataset sizes
- [ ] Different machines/environments
- [ ] Compared to previous versions

**Details:**


### 🧠 Memory Profiling
**If you've done memory profiling:**
```r
# Results from profvis::profvis() or similar tools
# Memory usage patterns you've observed
```

### ⚙️ Optimization Attempts
**What have you tried to improve performance?**
- [ ] Reduced dataset size
- [ ] Changed analysis parameters
- [ ] Increased available memory
- [ ] Used different methods
- [ ] Parallelization attempts
- [ ] Other: ___________

### 🔄 Workarounds
**Have you found any workarounds?**


### 📝 Additional Context
**Environment details:**
- Running on: [ ] Local machine [ ] Server [ ] Cloud platform [ ] HPC cluster
- Concurrent processes: 
- Available memory: 
- R session details: 

**Other context:**
Any other information that might help diagnose the performance issue.