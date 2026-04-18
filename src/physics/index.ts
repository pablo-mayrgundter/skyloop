import { RapierRigidBody, RigidBodyProps } from '@react-three/rapier'
import { useMemo } from 'react'

export type PhysicsBodyConfig = RigidBodyProps & {
  initialVelocity?: [number, number, number]
}

export const createGliderBody = (config: PhysicsBodyConfig) => ({
  type: 'dynamic' as const,
  linearVelocity: config.initialVelocity,
  ...config
})

export const createSectionBody = (config: PhysicsBodyConfig) => ({
  type: 'kinematicPosition' as const,
  ...config
})

export const stepPhysicsWorld = (deltaTime: number) => {
  // TODO: integrate rapier world stepping with global store
  // For now this is a placeholder to demonstrate API surface
  return deltaTime
}

export const usePhysicsBody = (config: PhysicsBodyConfig) => {
  const memoized = useMemo(() => config, [config])
  return memoized
}

export type PhysicsHandle = RapierRigidBody | null
