"use client";

import { Suspense, useCallback, useState } from "react";
import { Canvas } from "@react-three/fiber";
import { OrbitControls, ContactShadows } from "@react-three/drei";
import { HouseModel } from "./HouseModel";

interface HouseSceneProps {
  onPartSelect: (id: string | null) => void;
  onPartHover:  (id: string | null) => void;
}

export default function HouseScene({ onPartSelect, onPartHover }: HouseSceneProps) {
  const [selected, setSelected] = useState<string | null>(null);
  const [hovered, setHovered]   = useState<string | null>(null);

  const handleSelect = useCallback(
    (id: string | null) => { setSelected(id); onPartSelect(id); },
    [onPartSelect]
  );
  const handleHover = useCallback(
    (id: string | null) => { setHovered(id); onPartHover(id); },
    [onPartHover]
  );

  return (
    <Canvas
      shadows
      camera={{ position: [7, 5, 9], fov: 42 }}
      gl={{ antialias: true }}
      onPointerMissed={() => handleSelect(null)}
    >
      <color attach="background" args={["#F0FBF4"]} />

      {/* Lighting */}
      <ambientLight intensity={0.55} />
      <directionalLight
        position={[10, 12, 6]}
        intensity={1.3}
        castShadow
        shadow-mapSize={[2048, 2048]}
        shadow-camera-near={0.5}
        shadow-camera-far={60}
        shadow-camera-left={-12}
        shadow-camera-right={12}
        shadow-camera-top={12}
        shadow-camera-bottom={-12}
      />
      <directionalLight position={[-6, 4, -4]} intensity={0.3} />
      <hemisphereLight args={["#b9f0d0", "#5a8a5a", 0.4]} />

      <Suspense fallback={null}>
        <HouseModel
          selected={selected}
          hovered={hovered}
          onSelect={handleSelect}
          onHover={handleHover}
        />
        <ContactShadows
          position={[0, -0.4, 0]}
          opacity={0.35}
          scale={14}
          blur={2.5}
          far={10}
        />
      </Suspense>

      <OrbitControls
        target={[0, 1.5, 0]}
        minDistance={4}
        maxDistance={20}
        maxPolarAngle={Math.PI / 2.05}
        enablePan={false}
        enableDamping
        dampingFactor={0.06}
      />
    </Canvas>
  );
}
