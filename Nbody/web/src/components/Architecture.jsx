import React, { useEffect, useRef, useState } from 'react';

const Architecture = () => {
  const [isVisible, setIsVisible] = useState(false);
  const sectionRef = useRef(null);

  useEffect(() => {
    const observer = new IntersectionObserver(([entry]) => {
      if (entry.isIntersecting) {
        setIsVisible(true);
        observer.disconnect();
      }
    }, { threshold: 0.15 });
    
    if (sectionRef.current) {
      observer.observe(sectionRef.current);
    }
    return () => observer.disconnect();
  }, []);

  return (
    <section ref={sectionRef} id="architecture" className={`glass-section ${isVisible ? 'visible' : ''}`}>
      <h2>Technical Architecture</h2>
      <div className="card-grid three-col">
        <div className="card tech-card">
          <div className="icon">⚙️</div>
          <h3>Rust WGPU Engine</h3>
          <p>High-performance compute shaders executing on the Vulkan backend. Implements Velocity Verlet integration with batch command submission to prevent TDR timeouts during intense <code className="inline-math">O(N&sup2;)</code> gravity passes.</p>
        </div>
        <div className="card tech-card">
          <div className="icon">🔥</div>
          <h3>PyTorch ROCm</h3>
          <p>Vectorized Python implementations designed for AMD GPU tensor acceleration. Rapid prototyping for higher dimensions (up to 100D) using the <code>astro_env</code> Anaconda environment.</p>
        </div>
        <div className="card tech-card">
          <div className="icon">🐘</div>
          <h3>PostgreSQL Analytics</h3>
          <p>All snapshot geometries, phase spaces, and virial parameters are logged to a local PostgreSQL <code>hypercluster</code> database, easily queried via Julia and Python plotting suites.</p>
        </div>
      </div>
    </section>
  );
};

export default Architecture;
