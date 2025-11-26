import { Line, Sphere } from '@react-three/drei'
import { useMemo } from 'react'
import { Vector3 } from 'three'
import { EARTH_RENDER_RADIUS, LOS_ANGELES, SAN_FRANCISCO, interpolateGreatCircle } from '../config/simulation'

const Globe = ({ overlays }: { overlays: Record<string, boolean> }) => {
  const arcPoints = useMemo(() => {
    const samples: Vector3[] = []
    const steps = 64
    for (let i = 0; i <= steps; i += 1) {
      const t = i / steps
      samples.push(interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, t, EARTH_RENDER_RADIUS + 0.001))
    }
    return samples
  }, [])

  return (
    <group>
      <Sphere args={[EARTH_RENDER_RADIUS, 64, 64]}>
        <meshStandardMaterial
          color="#1d2f4f"
          roughness={0.8}
          metalness={0}
          emissive={overlays.showDayNight ? '#0c101b' : '#090c16'}
          emissiveIntensity={overlays.showDayNight ? 0.25 : 0.1}
        />
      </Sphere>
      <Line points={arcPoints} color="#7ad7f0" lineWidth={2} dashed dashSize={0.02} gapSize={0.02} />
      {overlays.showAtmosphere && (
        <Sphere args={[EARTH_RENDER_RADIUS * 1.05, 64, 64]}>
          <meshStandardMaterial color="#7bdcf5" transparent opacity={0.08} depthWrite={false} />
        </Sphere>
      )}
    </group>
  )
}

export default Globe
