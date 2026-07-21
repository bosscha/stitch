import React, { useEffect, useRef, useState } from 'react';

const PhysicsSection = () => {
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
    <section ref={sectionRef} id="physics" className={`glass-section ${isVisible ? 'visible' : ''}`}>
      <h2>The Physics of N-Dimensions</h2>
      <div className="card-grid">
        <div className="card">
          <h3>Generalized Gauss's Law</h3>
          <p>In standard 3D space, gravitational flux spreads across the surface of a sphere (<code className="inline-math">4&pi;R²</code>), resulting in the famous <code className="inline-math">1/R²</code> inverse-square law.</p>
          <p>In this project, we generalize Gauss's Law to <code className="inline-math">N</code>-dimensional space, where the surface area of a hypersphere scales proportionally to <code className="inline-math">R<sup>N-1</sup></code>. Consequently, our N-Body GPU shaders dynamically scale gravitational force mathematically as:</p>
          <div className="equation">
            <span className="math">F &propto; 1 / R<sup>N-1</sup></span>
          </div>
        </div>
        <div className="card feature-list">
          <h3>Dimensional Behaviors</h3>
          <ul>
            <li><span className="dim-badge">1D</span> <strong>Constant Force:</strong> <code className="inline-math">1/R&deg;</code> - Gravity does not decay with distance.</li>
            <li><span className="dim-badge">2D</span> <strong>Inverse Law:</strong> <code className="inline-math">1/R&sup1;</code> - Force scales logarithmically with potential.</li>
            <li><span className="dim-badge">3D</span> <strong>Inverse-Square:</strong> <code className="inline-math">1/R&sup2;</code> - Standard Newtonian gravity.</li>
            <li><span className="dim-badge">4D+</span> <strong>Hyperspatial:</strong> <code className="inline-math">1/R&sup3;</code> and beyond - Extreme short-range decay.</li>
          </ul>
        </div>
      </div>
    </section>
  );
};

export default PhysicsSection;
