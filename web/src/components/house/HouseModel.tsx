"use client";

import { useMemo } from "react";
import * as THREE from "three";
import type { ThreeEvent } from "@react-three/fiber";

export interface HousePartProps {
  selected: string | null;
  hovered: string | null;
  onSelect: (id: string | null) => void;
  onHover: (id: string | null) => void;
}

const PART_COLORS: Record<string, { base: string; hover: string; active: string }> = {
  roof:       { base: "#6B3A2A", hover: "#8B4D3A", active: "#D97706" },
  walls:      { base: "#D6CFC4", hover: "#C4BDB3", active: "#F59E0B" },
  windows:    { base: "#BAE6FD", hover: "#93C5FD", active: "#3B82F6" },
  door:       { base: "#1E293B", hover: "#2D3F55", active: "#818CF8" },
  foundation: { base: "#94A3B8", hover: "#7A90A8", active: "#10B981" },
  solar:      { base: "#1E3A5F", hover: "#2A4D7F", active: "#6366F1" },
  heatpump:   { base: "#CBD5E1", hover: "#B8C5D3", active: "#0EA5E9" },
};

function partColor(id: string, selected: string | null, hovered: string | null) {
  const c = PART_COLORS[id] ?? PART_COLORS.walls;
  if (selected === id) return c.active;
  if (hovered === id) return c.hover;
  return c.base;
}

function partEmissive(id: string, selected: string | null) {
  return selected === id ? (PART_COLORS[id]?.active ?? "#000000") : "#000000";
}

function handlers(
  id: string,
  selected: string | null,
  onSelect: (id: string | null) => void,
  onHover: (id: string | null) => void
) {
  return {
    onClick:      (e: ThreeEvent<MouseEvent>)  => { e.stopPropagation(); onSelect(selected === id ? null : id); },
    onPointerOver:(e: ThreeEvent<PointerEvent>) => { e.stopPropagation(); onHover(id); },
    onPointerOut: (e: ThreeEvent<PointerEvent>) => { e.stopPropagation(); onHover(null); },
  };
}

// ── Foundation ───────────────────────────────────────────────────────────────

function Foundation({ selected, hovered, onSelect, onHover }: HousePartProps) {
  const col = partColor("foundation", selected, hovered);
  const emit = partEmissive("foundation", selected);
  const h = handlers("foundation", selected, onSelect, onHover);
  return (
    <mesh position={[0, -0.2, 0]} castShadow receiveShadow {...h}>
      <boxGeometry args={[4.6, 0.4, 3.6]} />
      <meshStandardMaterial color={col} emissive={emit} emissiveIntensity={0.15} roughness={0.9} />
    </mesh>
  );
}

// ── Main walls ────────────────────────────────────────────────────────────────

function Walls({ selected, hovered, onSelect, onHover }: HousePartProps) {
  const col = partColor("walls", selected, hovered);
  const emit = partEmissive("walls", selected);
  const h = handlers("walls", selected, onSelect, onHover);
  return (
    <mesh position={[0, 1.25, 0]} castShadow receiveShadow {...h}>
      <boxGeometry args={[4, 2.5, 3]} />
      <meshStandardMaterial color={col} emissive={emit} emissiveIntensity={0.1} roughness={0.7} />
    </mesh>
  );
}

// ── Roof (triangular prism) ───────────────────────────────────────────────────

function Roof({ selected, hovered, onSelect, onHover }: HousePartProps) {
  const col = partColor("roof", selected, hovered);
  const emit = partEmissive("roof", selected);
  const h = handlers("roof", selected, onSelect, onHover);

  const shape = useMemo(() => {
    const s = new THREE.Shape();
    s.moveTo(-2.4, 0);
    s.lineTo(0, 1.5);
    s.lineTo(2.4, 0);
    s.lineTo(-2.4, 0);
    return s;
  }, []);

  const extOpts = useMemo(
    (): THREE.ExtrudeGeometryOptions => ({ depth: 3.6, bevelEnabled: false }),
    []
  );

  return (
    <mesh position={[0, 2.5, -1.8]} castShadow {...h}>
      <extrudeGeometry args={[shape, extOpts]} />
      <meshStandardMaterial
        color={col} emissive={emit} emissiveIntensity={0.15}
        roughness={0.8} side={THREE.DoubleSide}
      />
    </mesh>
  );
}

// ── Solar panels on right roof slope ─────────────────────────────────────────
// Right slope: from world (2.4, 2.5) to (0, 4.0) in XY.
// Slope angle from horizontal: atan2(1.5, 2.4) ≈ 32°.
// Flat slab normal after rotation.z = -atan2(1.5,2.4) aligns with slope outward normal.

function SolarPanels({ selected, hovered, onSelect, onHover }: HousePartProps) {
  const col = partColor("solar", selected, hovered);
  const emit = partEmissive("solar", selected);
  const h = handlers("solar", selected, onSelect, onHover);
  const tilt = -Math.atan2(1.5, 2.4); // match right-slope angle

  // 2×2 panel grid on the right slope
  const positions: [number, number, number][] = [
    [1.64, 3.54, -0.45],
    [1.64, 3.54,  0.45],
    [0.79, 3.01, -0.45],
    [0.79, 3.01,  0.45],
  ];

  return (
    <group {...h}>
      {positions.map((pos, i) => (
        <mesh key={i} position={pos} rotation={[0, 0, tilt]} castShadow>
          <boxGeometry args={[0.9, 0.05, 0.78]} />
          <meshStandardMaterial
            color={col} emissive={emit} emissiveIntensity={0.15}
            roughness={0.15} metalness={0.4}
          />
        </mesh>
      ))}
      {/* Panel grid lines */}
      {positions.map((pos, i) => (
        <mesh key={`g${i}`} position={[pos[0], pos[1] + 0.03, pos[2]]} rotation={[0, 0, tilt]}>
          <boxGeometry args={[0.88, 0.01, 0.01]} />
          <meshStandardMaterial color="#6366F1" />
        </mesh>
      ))}
    </group>
  );
}

// ── Window (reusable) ─────────────────────────────────────────────────────────

interface WindowProps extends HousePartProps {
  x: number; y: number; z: number; ry?: number;
}

function Window({ x, y, z, ry = 0, selected, hovered, onSelect, onHover }: WindowProps) {
  const h = handlers("windows", selected, onSelect, onHover);
  const active = selected === "windows";
  const isHovered = hovered === "windows";
  const frameCol = active ? "#3B82F6" : isHovered ? "#93C5FD" : "#E5E7EB";
  const glassCol = active ? "#BFDBFE" : "#BAE6FD";
  return (
    <group position={[x, y, z]} rotation={[0, ry, 0]} {...h}>
      {/* Frame */}
      <mesh>
        <boxGeometry args={[0.84, 0.74, 0.07]} />
        <meshStandardMaterial color={frameCol} roughness={0.3} metalness={0.2} />
      </mesh>
      {/* Glass */}
      <mesh position={[0, 0, 0.04]}>
        <boxGeometry args={[0.70, 0.60, 0.04]} />
        <meshStandardMaterial color={glassCol} roughness={0.05} metalness={0.1} transparent opacity={0.75} />
      </mesh>
      {/* Vertical mullion */}
      <mesh position={[0, 0, 0.075]}>
        <boxGeometry args={[0.04, 0.60, 0.025]} />
        <meshStandardMaterial color={frameCol} />
      </mesh>
      {/* Horizontal rail */}
      <mesh position={[0, 0, 0.075]}>
        <boxGeometry args={[0.70, 0.04, 0.025]} />
        <meshStandardMaterial color={frameCol} />
      </mesh>
    </group>
  );
}

// ── Front door ────────────────────────────────────────────────────────────────

function Door({ selected, hovered, onSelect, onHover }: HousePartProps) {
  const col = partColor("door", selected, hovered);
  const emit = partEmissive("door", selected);
  const h = handlers("door", selected, onSelect, onHover);
  const frameCol = selected === "door" ? "#818CF8" : "#CBD5E1";
  return (
    <group>
      {/* Frame */}
      <mesh position={[0, 1.0, 1.52]} {...h}>
        <boxGeometry args={[1.12, 2.12, 0.05]} />
        <meshStandardMaterial color={frameCol} roughness={0.5} />
      </mesh>
      {/* Door leaf */}
      <mesh position={[0.02, 0.98, 1.57]} castShadow {...h}>
        <boxGeometry args={[0.94, 1.94, 0.08]} />
        <meshStandardMaterial color={col} emissive={emit} emissiveIntensity={0.12} roughness={0.4} />
      </mesh>
      {/* Small window in door */}
      <mesh position={[0.02, 1.5, 1.63]}>
        <boxGeometry args={[0.5, 0.4, 0.04]} />
        <meshStandardMaterial color="#BAE6FD" transparent opacity={0.7} roughness={0.05} />
      </mesh>
      {/* Handle */}
      <mesh position={[-0.28, 0.98, 1.66]}>
        <boxGeometry args={[0.04, 0.18, 0.04]} />
        <meshStandardMaterial color="#B8A060" metalness={0.8} roughness={0.2} />
      </mesh>
      {/* Step */}
      <mesh position={[0, -0.1, 1.86]} receiveShadow>
        <boxGeometry args={[1.4, 0.2, 0.5]} />
        <meshStandardMaterial color="#94A3B8" roughness={0.9} />
      </mesh>
    </group>
  );
}

// ── Heat pump unit ────────────────────────────────────────────────────────────

function HeatPump({ selected, hovered, onSelect, onHover }: HousePartProps) {
  const col = partColor("heatpump", selected, hovered);
  const emit = partEmissive("heatpump", selected);
  const h = handlers("heatpump", selected, onSelect, onHover);
  return (
    <group position={[2.55, 0.45, 0.6]} {...h}>
      <mesh castShadow>
        <boxGeometry args={[0.3, 0.9, 0.65]} />
        <meshStandardMaterial color={col} emissive={emit} emissiveIntensity={0.15} roughness={0.4} metalness={0.2} />
      </mesh>
      {/* Grille slats */}
      {[-0.15, 0, 0.15].map((dz, i) => (
        <mesh key={i} position={[0.16, 0.05, dz]}>
          <boxGeometry args={[0.02, 0.6, 0.12]} />
          <meshStandardMaterial color="#9CA3AF" roughness={0.6} />
        </mesh>
      ))}
      {/* Base */}
      <mesh position={[0, -0.52, 0]}>
        <boxGeometry args={[0.36, 0.14, 0.72]} />
        <meshStandardMaterial color="#6B7280" roughness={0.8} />
      </mesh>
    </group>
  );
}

// ── Chimney (decorative) ──────────────────────────────────────────────────────

function Chimney() {
  return (
    <group>
      <mesh position={[-0.7, 3.9, -0.5]} castShadow>
        <boxGeometry args={[0.44, 1.4, 0.44]} />
        <meshStandardMaterial color="#64748B" roughness={0.9} />
      </mesh>
      <mesh position={[-0.7, 4.66, -0.5]}>
        <boxGeometry args={[0.54, 0.12, 0.54]} />
        <meshStandardMaterial color="#475569" roughness={0.8} />
      </mesh>
    </group>
  );
}

// ── Ground ────────────────────────────────────────────────────────────────────

function Ground({ onSelect }: { onSelect: (id: string | null) => void }) {
  return (
    <mesh
      position={[0, -0.41, 0]}
      rotation={[-Math.PI / 2, 0, 0]}
      receiveShadow
      onClick={() => onSelect(null)}
    >
      <planeGeometry args={[40, 40]} />
      <meshStandardMaterial color="#4A7A54" roughness={1} />
    </mesh>
  );
}

// ── Assembled house ───────────────────────────────────────────────────────────

export function HouseModel({ selected, hovered, onSelect, onHover }: HousePartProps) {
  const p = { selected, hovered, onSelect, onHover };
  return (
    <group>
      <Ground onSelect={onSelect} />
      <Foundation {...p} />
      <Walls {...p} />
      <Roof {...p} />
      <SolarPanels {...p} />
      <Chimney />

      {/* Front windows */}
      <Window x={-1.1} y={1.55} z={1.52} {...p} />
      <Window x={1.35} y={1.55} z={1.52} {...p} />

      {/* Side windows */}
      <Window x={-2.02} y={1.55} z={0.0}  ry={-Math.PI / 2} {...p} />
      <Window x={ 2.02} y={1.55} z={-0.4} ry={ Math.PI / 2} {...p} />

      {/* Back windows */}
      <Window x={-0.9} y={1.55} z={-1.52} ry={Math.PI} {...p} />
      <Window x={ 0.9} y={1.55} z={-1.52} ry={Math.PI} {...p} />

      <Door {...p} />
      <HeatPump {...p} />
    </group>
  );
}
