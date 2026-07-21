import React from 'react';
import Hero from './components/Hero';
import PhysicsSection from './components/PhysicsSection';
import IncomingObject from './components/IncomingObject';
import Architecture from './components/Architecture';

function App() {
  return (
    <>
      <Hero />
      <main>
        <PhysicsSection />
        <IncomingObject />
        <Architecture />
      </main>
      <footer>
        <p>&copy; 2026 GAIA Hypercluster Research | Simulated locally via Antigravity</p>
      </footer>
    </>
  );
}

export default App;
