import js from '@eslint/js'
import tseslint from 'typescript-eslint'
import react from 'eslint-plugin-react'
import reactHooks from 'eslint-plugin-react-hooks'
import globals from 'globals'

export default tseslint.config(
  js.configs.recommended,
  ...tseslint.configs.recommended,
  {
    files: ['**/*.{js,jsx,ts,tsx}'],
    languageOptions: {
      globals: {
        ...globals.browser,
        ...globals.node,
        ...globals.es2021
      }
    },
    plugins: {
      react: react,
      'react-hooks': reactHooks
    },
    settings: {
      react: {
        version: 'detect'
      }
    },
    rules: {
      ...react.configs.recommended.rules,
      ...reactHooks.configs.recommended.rules,
      'semi': ['error', 'never'],
      'quote-props': ['error', 'as-needed'],
      'indent': ['error', 2],
      'max-len': ['error', 120],
      'object-curly-spacing': ['error', 'always'],
      'react/react-in-jsx-scope': 'off',
      'react/no-unknown-property': ['error', { ignore: ['args', 'attach', 'intensity', 'position', 'castShadow', 'emissive', 'wireframe', 'object'] }],
      '@typescript-eslint/no-unused-vars': ['error', { argsIgnorePattern: '^_' }]
    }
  },
  {
    ignores: ['node_modules/**', 'dist/**', 'build/**', '*.config.js', '*.config.ts']
  }
)

